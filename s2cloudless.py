import ee


S2Cloudless = ee.ImageCollection

class S2CloudProbabilityCollection(ee.ImageCollection):
    def __init__(self) -> None:
        """a class to represent the S2CloudProbability Image Collection
        Extends: ee.ImageCollection
        """
        self._arg = "COPERNICUS/S2_CLOUD_PROBABILITY"
        super().__init__(self._arg)

    @property
    def arg(self):
        return self._arg


class S2SRCollection(ee.ImageCollection):
    """ a class to represent the S2SR Image Collection"""
    def __init__(self):
        """ Extends: ee.ImageCollection """
        self._arg = "COPERNICUS/S2_SR"
        super().__init__(self._arg)
        
    @property
    def arg(self):
        return self._arg


class S2CloudCollection(ee.ImageCollection):
    def __init__(self, cld_prb: S2CloudProbabilityCollection, sr_sr: S2SRCollection) -> None:
        """ Extends: ee.ImageCollection """
        self.cld_prb = cld_prb
        self.sr_sr = sr_sr
        super().__init__(self._join(self.sr_sr, self.cld_prb))
    
    def _join(self, primary, secondary)-> ee.ComputedObject:
        """joing two image collections together on the system:index property

        Args:
            primary (S2SRCollection): s2_sr collection
            secondary (S2CloudProbabilityCollection): cloud probability collection

        Returns:
            ee.ComputedObject: the joined collection
        """
        return ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': primary,
        'secondary': secondary,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
            })
        })


class S2CloudlessBuilder:
    def __init__(self) -> None:
        self._col = None
    
    @property
    def product(self) -> S2Cloudless:
        return self._col

    @product.setter
    def product(self, col: S2Cloudless) -> None:
        self._col = col
    
    def add_cloud_bands(self, cld_prb_thresh: int = 60) -> None:
        def wrapper(img):
            # Get s2cloudless image, subset the probability band.
            cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

            # Condition s2cloudless by the probability threshold value.
            is_cloud = cld_prb.gt(cld_prb_thresh).rename('clouds')

            # Add the cloud probability layer and cloud mask as image bands.
            return img.addBands(ee.Image([cld_prb, is_cloud]))
        self.col = self.col.map(wrapper)
        return self
    
    def add_shadow_bands(self, nir_drk_thresh: float = 0.15, cld_prj_dist: int = 2) -> None:
        def wrapper(img):
            # Identify water pixels from the SCL band.
            not_water = img.select('SCL').neq(6)

            # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
            SR_BAND_SCALE = 1e4
            dark_pixels = img.select('B8').lt(nir_drk_thresh*SR_BAND_SCALE)\
                .multiply(not_water).rename('dark_pixels')

            # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
            shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))

            # Project shadows from clouds for the distance specified by the cld_prj_dist input.
            cld_proj = (img.select('clouds')\
                .directionalDistanceTransform(shadow_azimuth, cld_prj_dist*10)
                .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
                .select('distance')
                .mask()
                .rename('cloud_transform'))

            # Identify the intersection of dark pixels with cloud shadow projection.
            shadows = cld_proj.multiply(dark_pixels).rename('shadows')

            # Add dark pixels, cloud projection, and identified shadows as image bands.
            return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))
        self.col = self.col.map(wrapper)
        return self

    def add_cld_shdw_mask(self, buffer: int = 100) -> None:
        def wrapper(img):
            # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
            is_cld_shdw = img.select('clouds').add(img.select('shadows')).gt(0)

            # Remove small cloud-shadow patches and dilate remaining pixels by buffer input.
            # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
            is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(buffer*2/20)
                .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
                .rename('cloudmask'))

            # Add the final cloud-shadow mask to the image.
            return img.addBands(is_cld_shdw)
        self.col = self.col.map(wrapper)
        return self

    def apply_cld_shdw_mask(self):
        def wrapper(img):
            # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
            not_cld_shdw = img.select('cloudmask').Not()

            # Subset reflectance bands and update their masks, return the result.
            return img.select('B.*').updateMask(not_cld_shdw)
        self.col = self.col.map(wrapper)
        return self

    def build(self) -> S2Cloudless:
        return self.col
    

class S2CloudlessDirector:
    def __init__(self) -> None:
        self._builder = None
    
    @property
    def builder(self) -> S2CloudlessBuilder:
        return self._builder

    @builder.setter
    def builder(self, builder: S2CloudlessBuilder) -> None:
        self._builder = builder
    
    def build(self, cld_prb_thresh: int = 50, nir_drk_thresh: float = 0.15, cld_prj_dist =1,
              buffer: int = 50) -> S2Cloudless:
        self._builder\
            .add_cloud_bands(cloud_prob_thresh=cld_prb_thresh)\
            .add_shadow_bands()\
            .add_cld_shdw_mask()\
            .apply_cld_shdw_mask()\
            .build()


def build_s2_cloudless(aoi: ee.Geometry, date_range: tuple[str], cloud_filter: int = 60, 
                       cld_prb_thresh: int = 60, nir_drk_thresh: float = 0.15, cld_prj_dst = 2,
                       buffer = 100) -> S2Cloudless:
    """ executes the build work flow for the S2Cloudless class
    
    inspired by https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless
    
    Parameters:
        aoi (ee.Geometry): area of interest
        date_range (tuple[str]): start and end date in the format 'YYYY-MM-DD'
        cloud_filter (int): maximum cloud cover percentage
        cld_prb_thresh (int): cloud probability threshold
        nir_drk_thresh (float): dark NIR threshold
        cl_prj_dst (int): cloud projection distance
        buffer (int): buffer distance
    Returns:
        S2Cloudless: S2Cloudless object
    """
    s2_sr = S2SRCollection()\
        .filterBounds(aoi)\
        .filterDate(*date_range)\
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_filter))

    s2_cld_prb = S2CloudProbabilityCollection()\
        .filterBounds(aoi)\
        .filterDate(*date_range)
    
    s2_cld_col = S2CloudCollection(s2_cld_prb, s2_sr)
    
    bldr = S2CloudlessBuilder()
    bldr.col = s2_cld_col
    
    dctr = S2CloudlessDirector()
    dctr.builder = bldr
    dctr.build(cld_prb_thresh, nir_drk_thresh, cld_prj_dst, buffer)
    
    return bldr.col
