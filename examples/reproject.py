import json
import sys
import math

# Opening JSON file
f = open(sys.argv[1])

# returns JSON object as
# a dictionary
data = json.load(f)

def reproject(point):
	IRAD = 180.0 / math.pi

	x = point[0]
	y = point[1]
	lat = (1.5707963267948966 - (2.0 * math.atan(math.exp(-y / 6378137.0)))) * IRAD
	lon = x / 111319.4907932735677
	return [lon, lat]

# Iterating through the json
# list
for i in data['features']:
	if i['geometry']['type'] == 'Point':
		i['geometry']['coordinates'] = reproject(i['geometry']['coordinates'])
		#  print(i)
	if i['geometry']['type'] == 'LineString':
		for c in range(len(i['geometry']['coordinates'])):
			i['geometry']['coordinates'][c] = reproject(i['geometry']['coordinates'][c])

# Closing file
f.close()

# Serializing json
json_object = json.dumps(data, indent=2)

# Writing to sample.json
with open(sys.argv[1] + ".new", "w") as outfile:
    outfile.write(json_object)
