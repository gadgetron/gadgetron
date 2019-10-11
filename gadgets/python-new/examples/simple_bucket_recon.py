import gadgetron
from gadgetron.util import cifftn,cfftn
import time
import numpy as np
import ismrmrd
from itertools import starmap
import logging


def build_buffers(buckets):

    for id, bucket in buckets:
        assert(type(bucket) == gadgetron.AcquisitionBucket)
        buffer = np.stack([acq.data for acq in bucket.data],axis=2)
        yield bucket.data[0],buffer

def reconstruct_image(acquisition,buffer):
    return acquisition, gadgetron.util.cifftn(buffer, axes=[1,2])

def combine_channels(acq, image):
    return acq, np.sqrt(np.sum(np.abs(image)**2,axis=0))

def create_ismrmrd_images(acquisition, image):
    return ismrmrd.image.Image.from_array(image, acquisition=acquisition)

@gadgetron.external.module
def bucket_recon(connection):

    logging.info("Python reconstruction running.")
    start = time.time()
    connection.filter(gadgetron.legacy.AcquisitionBucket)
    buffers = build_buffers(connection)

    images = starmap(reconstruct_image,buffers)
    images = starmap(combine_channels,images)
    images = starmap(create_ismrmrd_images,images)

    for image in images:
        connection.send(image)

    logging.info(f"Python reconstruction done. Duration: {(time.time() - start):.2f} s")

def main():
    import gadgetron.external
    logging.basicConfig(level=logging.INFO)
    gadgetron.external.listen(9100, bucket_recon)


if __name__ == '__main__':
    main()
