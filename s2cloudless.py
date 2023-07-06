import ee

S2CloudCollection = ee.ImageCollection
S2CloudProbability = ee.ImageCollection
S2Cloudless = ee.ImageCollection
S2Collection = ee.ImageCollection



class S2CloudlessBuilder:
    def __init__(self) -> None:
        self.col = None
    
    def add_cloud_bands(self, cloud_prob_thresh: int = 60) -> None:
        def wrapper(img):
            # Get s2cloudless image, subset the probability band.
            cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

            # Condition s2cloudless by the probability threshold value.
            is_cloud = cld_prb.gt(cloud_prob_thresh).rename('clouds')

            # Add the cloud probability layer and cloud mask as image bands.
            return img.addBands(ee.Image([cld_prb, is_cloud]))
        self.col = self.col.map(wrapper)
        return self
    
    def add_shadow_bands(self, NIR_DRK_THRESH: float = 0.15, CLD_PRJ_DIST: int = 2) -> None:
        def wrapper(img):
            # Identify water pixels from the SCL band.
            not_water = img.select('SCL').neq(6)

            # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
            SR_BAND_SCALE = 1e4
            dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')

            # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
            shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

            # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
            cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
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

    def add_cld_shdw_mask(self, BUFFER: int = 100) -> None:
        def wrapper(img):
            # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
            is_cld_shdw = img.select('clouds').add(img.select('shadows')).gt(0)

            # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
            # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
            is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
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
    


def build_s2_cloudless(aoi: ee.Geometry, date_range: tuple[str], cloud_filter: int = 60, 
                       cld_prb_thresh: int = 60, nir_drk_thresh: float = 0.15, cld_prj_dst = 2,
                       buffer = 100) -> S2Collection:
    """ Builds a S2Cloudless object. """
    s2_cld_col = S2CloudCollection('COPERNICUS/S2_SR')\
            .filterBounds(aoi)\
            .filterDate(*date_range)\
            .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', cloud_prob))
        
    s2_cld_prob = S2CloudProbability('COPERNICUS/S2_CLOUD_PROBABILITY')\
        .filterBounds(aoi)\
        .filterDate(*date_range)\
        
    cld_prb = S2Cloudless(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': s2_cld_col,
        'secondary': s2_cld_prob,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))
    
    bldr = S2CloudlessBuilder()
    bldr.col = cld_prb
    
    bldr.add_cloud_bands()
    bldr.add_shadow_bands()
    bldr.add_cld_shdw_mask()
    bldr.apply_cld_shdw_mask()
    bldr.build()
    
    return bldr.col
