SpaceNet Buildings Exploration

Transform SpaceNet geojson buidling labels data into raster masks.
Download data via:
    aws s3api get-object --bucket spacenet-dataset \
    --key AOI_1_Rio/processedData/processedBuildingLabels.tar.gz \
    --request-payer requester processedBuildingLabels.tar.gz

Download spacenet utilities from:
   https://github.com/SpaceNetChallenge/utilities/tree/master/python/spaceNet 

For further details, see:
    https://medium.com/the-downlinq/getting-started-with-spacenet-data-827fd2ec9f53
