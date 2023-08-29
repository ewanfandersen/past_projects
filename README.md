# This pipeline allows its user to approximate the age of a galaxy cluster using theoretical isochrones and Gaia Observatory data from the VizeiR database. 
# The isochrones in question are theoretical models of the shape that a galaxy cluster will make when plotted on a Color-Magnitude diagram, of which this code will take five.
# First, the target galaxy cluster is isolated with filters excluding data points outside a set of parameters.
# Next, the isolated cluster data is plotted onto a Color-Magnitude diagram.
# Finally, the five isochrones are plotted onto the same Color-Magnitude diagram and the age can be approximated.
# The pipeline is set up for the galaxy cluster NGC 2169 as an example
