#!/usr/bin/python

import os; 
import sys; 
import math;
import pdb;
import time;
import string;
import imp;
pygeo = imp.load_source('pygeocoder', 'pygeocoder/pygeocoder.py')


from datetime import  datetime, tzinfo, timedelta; 

class UTC(tzinfo):
	def utcoffset(self, dt):
		return timedelta(0)
	def tzname(self,dt):
		return "UTC"
	def dst(self, dt):
		return timedelta(0)


class EST(tzinfo):
	def utcoffset(self, dt):
		return timedelta(hours=-5)
	def tzname(self,dt):
		return "EST"
	def dst(self, dt):
		return timedelta(0)

#$GPGGA,043440.600,3904.1210,N,07656.6711,W,2,9,0.94,90.3,M,-33.5,M,0000,0000*55
#$GPRMC,043440.600,A,3904.1210,N,07656.6711,W,0.07,276.64,130710,,,D*75

M_PER_MI=1609.344 #meters/mi
KNOTS_PER_MPH=0.868976242

# helper functions
def to_radian (deg):
	return (math.pi*deg)/180.0

def compute_distance (lata,lona,latb,lonb,alta,altb,use_altitude=False):
	EARTH_RADIUS=6378137.0 #meters
	lata = to_radian(lata)
	latb = to_radian(latb)
	lonb = to_radian(lonb)
	lona = to_radian(lona)
	dlat = latb-lata
	dlon = lonb-lona;
	dalt = altb-alta;
#Haversine formula
	angle = math.sqrt((math.sin(dlat/2.0)**2) + (math.cos(lata)*math.cos(latb)*(math.sin(dlon/2)**2)))
	distance = (EARTH_RADIUS * 2.0 * math.asin(angle))
#pythagorean theorem to account for elevation
	if (use_altitude and dalt > 0.0):
		distance = math.sqrt(distance**2 + dalt**2)
#	distance /= M_PER_MI
	return distance

def minutes_to_degrees(coord, direction):
	coord_float = coord/100.0;
	coord_degree = math.floor(coord_float)+(coord_float%1.0)/0.6 # modulo 1 gives you just the fractional part, 0.6 = 60min/100
	if direction =="S" or direction == "W":
		coord_degree = coord_degree * -1.0;
	return coord_degree
def nmea_date_and_time_to_datetime(nmea_date, nmea_time):
	hour = int(nmea_time[0:2]); 
	minute = int(nmea_time[2:4]);
	second = int(nmea_time[4:6]);
	microsecond = int(nmea_time[7:10])*1000;

	day = int(nmea_date[0:2])
	month = int(nmea_date[2:4]);
	year = 2000 + int(nmea_date[4:6]);
	#XXX: this certainly would have been a lot easier with -timedelta(hours=5) :|  
	return datetime(year, month, day, hour, minute, second, microsecond, tzinfo=UTC()).astimezone(EST())

def valid_gps_checksum(sentence):
	split_sentence = sentence.split("*")
	if len(split_sentence) < 2 or len(split_sentence) > 12:
		return False; 
	checksum=0;
	for c in split_sentence[0][1:]:
		checksum = checksum ^ ord(c)
	return hex(checksum)[2:].upper() == split_sentence[1];

class GGASentence: 
	def __init__(self, sentence_fields):
		assert sentence_fields[0] == "$GPGGA", "Not a GGA Sentence: %s" % ",".join(sentence_fields)
		self.timestamp_string = sentence_fields[1]
		lat,latd = float(sentence_fields[2]),sentence_fields[3]
		lon,lond = float(sentence_fields[4]),sentence_fields[5] 
		pos_fix,num_sats = int(sentence_fields[6]),int(sentence_fields[7])
		altitude = float(sentence_fields[9])

		self.lat = minutes_to_degrees(lat, latd);
		self.lon = minutes_to_degrees(lon, lond);
		self.altitude = altitude
		self.pos_fix = pos_fix
		self.num_sats = num_sats
	def __repr__(self):
		return ("[%s] GGA: %f,%f alt=%f m, %d satellites (%d)" % (self.timestamp_string, self.lat, self.lon, self.altitude, self.num_sats, self.pos_fix))

class RMCSentence:
	def __init__(self, sentence_fields):
		assert sentence_fields[0] == "$GPRMC", "Not an RMC Sentence: %s" % ",".join(sentence_fields)
#		print("RMC: %s" % ",".join(sentence_fields))
		self.timestamp_string = sentence_fields[1]
		self.velocity=float(sentence_fields[7])/KNOTS_PER_MPH
#		print ("velocity is %f (%s)" % (self.velocity,sentence_fields[7]));
		self.bearing=float(sentence_fields[8])
		self.datestamp_string=sentence_fields[9]
	def __repr__(self):
		return ("[%s] RMC: v=%fm/s bearing=%f %s" % (self.timestamp_string, self.velocity, self.bearing, self.datestamp_string));

class GPSPoint:
	def __init__(self, gga, rmc): 
		self.lat = gga.lat
		self.lon = gga.lon
		self.timestamp = nmea_date_and_time_to_datetime(rmc.datestamp_string, rmc.timestamp_string);
		self.velocity = rmc.velocity
		self.altitude = gga.altitude
	def distance_to(self, other, use_miles=False):
		distance = compute_distance(self.lat,self.lon, other.lat, other.lon, self.altitude, other.altitude);
		if (distance == 0):
			pass
#			sys.stderr.write("Zero distance between %s and %s\n" % (self, other));
		if (use_miles):
			distance = distance / M_PER_MI;
		return distance
	def get_reverse_geocoding(self, index, spaces=False):
		#where index is the index into the result string: 2 is the city, 0 is the full address
		result = pygeo.Geocoder.reverse_geocode(self.lat, self.lon);
		print "----------"+str(result[index])
		return_string = str(result[index]).replace(", USA","")
		if not spaces:
			return_string = return_string.replace(" ", "_")
			return_string = return_string.replace(",", "" )
		return return_string;
		
	def __repr__(self):
		return "[%s] %f,%f v=%fmph alt=%fm" % (self.timestamp.strftime("%m.%d.%Y %H:%M"), self.lat, self.lon, self.velocity, self.altitude)
	def javascript_string(self):
		return "new google.maps.LatLng(%f, %f)" % (self.lat, self.lon);

class GPSSegment:
	def __init__(self, points_list):
		self.points = points_list;
	def length(self, use_miles=False):
		distance = 0.0;
		for index in range(len(self.points)-1):
			distance += self.points[index].distance_to(self.points[index+1]);
		if use_miles:
			distance /= M_PER_MI;
		return distance; 
	def time_taken(self, in_hours=False):
		delta = (self.points[-1].timestamp - self.points[0].timestamp).seconds;
		if in_hours:
			return delta/3600.0
		else:
			return delta

	def num_points(self):
		return len(self.points);
	def max_speed(self):
		max_speed = self.points[0].velocity;
		for p in self.points:
			max_speed = max(p.velocity, max_speed);
		return max_speed;
	def max_altitude(self):
		max_altitude = self.points[0].altitude;
		for p in self.points:
			max_altitude = max(p.altitude, max_altitude);
		return max_altitude;

	def min_altitude(self):
		min_altitude = self.points[0].altitude;
		for p in self.points:
			min_altitude = min(p.altitude, min_altitude);
		return min_altitude;

	def average_speeds(self):
			return {"moving": self.length(True)/self.time_taken(True), "overall":self.length(True)/self.time_taken(True)};
	def slice_velocity_range(self,low,high):
		ins=[];
		outs=[];

		in_segments = [];
		out_segments = [];
		starts = [];
		ends = [];
		last_one_in_segment=-1;
		for i,p in enumerate(self.points):
			#if this point doesn't fit
			if not (p.velocity > low and p.velocity <= high):
				#but if the last one did
				if (last_one_in_segment > -1):
					#this is the end of a segment
					ins[-1] = (last_one_in_segment, i)
#					ends.append(i)
				last_one_in_segment=-1;
			#but, if the point DOES fit AND the last one did not
			elif last_one_in_segment<0:
				# start up a new segment
				ins.append((i, -1))
#				starts.append(i);
				last_one_in_segment=i;
		if len(ins) == 0:
			return [], [self]; # everyone is out, no one is in

		first_start, first_stop = ins[0];
		last_start, last_stop = ins[-1];

		#if the first isn't zero, add an out segment
		if (first_start != 0):
			outs.append((0,first_start)); 

		# last one has a start but no end, so extend the last in to go to the end of the range
		if last_stop == -1:
			last_begin, last_end = ins[-1]	
			ins[-1] = (last_begin, len(self.points))

		# but if the last stop is less than the last point, then the segment is out
		elif last_stop < len(self.points):
			outs.append((last_stop,len(self.points)))
		# now handle all the normal cases
		for x in range(len(ins)-1):
			start,end = ins[x]
			next_start, next_end = ins[x+1]
			outs.append((end,next_start));
			
#		print "%s:%d\n%s:%d\n" % (ins, len(ins), outs, len(outs));
		for i in ins:
			begin,end = i
			in_segments.append(GPSSegment(self.points[begin:end]));
		for o in outs:
			begin,end = o
			out_segments.append(GPSSegment(self.points[begin:end]));

#			# the else condition here is just that the point does fit and so did the last one, so there is nothing to do but continue
#
##TODO: all of this is a terrible mess, need some much cleaner way to do this
#
#		# if we didn't start any segments, this whole segment is out
#		if len(starts) == 0:
#			return [], [self]
#		# if we did start a segment, but didn't end -- then this whole segment is in 
#		if len(ends) == 0:
#			return [self], []
#
#		if (ends[0] != 0):
#			sys.stderr.write("prepending zero\n");
#			ends = [0] + ends; #prepend a zero 
#
#		if (len(starts) < len(ends)):
#			sys.stderr.write("appending last point %d\n" % len(self.points));
#			starts.append(len(self.points));
#		elif (len(starts) == len(ends)):
#			sys.stderr.write("appending to ends %d\n" % len(self.points));
#			ends.append(len(self.points));
#			
#
#		sys.stderr.write ("%s:%d\n%s:%d \n" % (starts,len(starts),ends,len(ends)))
##		assert len(starts) == len(ends), "unequal starts, ends %s %s" %( starts, ends )
#		for i in range(len(ends)):
#			if (i < len(ends)):
#				out_segments.append(GPSSegment(self.points[ends[i]:starts[i]]))
#				sys.stderr.write ("i=%d out segment from %d to %d\n" % (i,ends[i],starts[i]))
#
#
#			# don't do this for the last iteration since it will be out of bounds
#			if (i < len(starts)-1):
#				in_segments.append(GPSSegment(self.points[starts[i]:ends[i+1]]))
#				sys.stderr.write ("i=%d: in segment from %d to %d\n" % (i,starts[i],ends[i+1]))
		combined_points = 0;
		for i in in_segments+out_segments:
			combined_points += i.num_points();
		assert combined_points == self.num_points(), "Lost some points (%d originally, now %d)" %(self.num_points(),combined_points)
		return in_segments, out_segments
		


	def slice_largest_constant_velocity_segment(self,vrange):
		segments = range(len(self.points))
		start,end = 0,0;
		
		for i,p in enumerate(self.points):
			segments[i] = 0
			start_range, end_range = p.velocity-vrange/2.0, p.velocity+vrange/2.0
			for j in range(i,len(self.points)):
				if self.points[j].velocity < end_range and self.points[j].velocity > start_range:
					segments[i] += 1;
				else:
					break;
		#print segments;
		max_index=0
		max_length=0;
		for x,length in enumerate(segments):
			if length > max_length:
				max_length = length;
				max_index = x;
		sys.stderr.write("Max segment is from %d to %d (len=%d/%d)\n" % (max_index, max_index+max_length, max_length, len(self.points)));
		return GPSSegment(self.points[0:max_index]), GPSSegment(self.points[max_index:max_index+max_length]), GPSSegment(self.points[max_index+max_length:])

				
	def javascript_string(self, line_color="#FF0000"):
		js_string = "{ \"points\" : ["
		for index,point in enumerate(self.points):
			if index < len(self.points)-1:
				js_string += "\t%s,\n" % point.javascript_string()
			else:
				js_string += "\t%s \n" % point.javascript_string()
	
		js_string += "], \"color\": \"%s\" }" % line_color
		return js_string;

	def summary_string(self, start_str="", end_str="\n"):
		if len(self.points) == 0:
				return "";
		start_point, stop_point = self.points[0],self.points[-1]
		average_speeds = self.average_speeds()
		str = "%sStarted=%s, Ended=%s%s" % (start_str, start_point.timestamp,stop_point.timestamp,end_str)
		str += "%sTrip Distance=%fm, %fmi, Trip Time=%s%s" % (start_str, self.length(),self.length(use_miles=True), stop_point.timestamp - start_point.timestamp,end_str)
		str += "%sMax Speed=%fmph, Moving Average Speed=%f, Overall Average=%f%s" % (start_str, self.max_speed(), average_speeds["moving"], average_speeds["overall"],end_str)
		str += "%sAltitude difference: %d - %d%s" % (start_str,self.min_altitude(), self.max_altitude(),end_str);
		return str;
			
	def __repr__(self):
		return self.summary_string;

		
if __name__ == "__main__":

	if len(sys.argv) < 2:
		print ("Not enough arguments");
	
	gps_file = open(sys.argv[1], "r"); 
	lines = []
	points = []
	total_points = 0; 
	added_points = 0; 
	while 1:
		next_line = gps_file.readline();
		
		if not next_line: 
			break;
		next_line=next_line.rstrip();
		next_line = next_line[string.rfind(next_line,"$"):]
		if not valid_gps_checksum(next_line):
			sys.stderr.write("Throwing away invalid checksum line '%s'\n" % next_line);
			continue;

		lines.append(next_line.split(","))
		if len(lines) < 2:
			continue;
		if lines[0][1] == lines[1][1] and lines[0][0] != lines[1][0]:
#			print ("Valid line '%s' and '%s'" % (",".join(lines[0]), ",".join(lines[1])));
			if (lines[0][0] == "$GPGGA"):
				gga = GGASentence(lines[0]);
				rmc = RMCSentence(lines[1]);
			else:
				rmc = RMCSentence(lines[0]);
				gga = GGASentence(lines[1]);

			if (gga.num_sats > 5):
				points.append (GPSPoint(gga,rmc));
				added_points = added_points+1;
			total_points = total_points+1;


			del(lines[0:1])
		else:
			del(lines[0])
			continue;
	starting_point = points[0].get_reverse_geocoding(2)
	ending_point = points[-1].get_reverse_geocoding(2)
	print str(starting_point +" to "+ ending_point)
	js_filename = "html/"
	js_filename += points[0].timestamp.strftime("%Y_%b_%d_%H.%M");
	js_filename += "_%s__%s.js"%(starting_point,ending_point)
	js_file = open(js_filename,"w"); 
	print "Saving javascriptto '"+js_filename+"'";


	sys.stderr.write("Original points=%d, added=%d\n" % (total_points,added_points));
	car_ranges = [35,55,55,100]
	bike_ranges = [5,12,12,40]
	run_ranges = [3,5,5,10]

	if len(sys.argv) == 3 and sys.argv[2] == "bike":
		ranges = bike_ranges
	else:
		ranges = car_ranges

	segment = GPSSegment(points);	
	fast_segments,out_segments = segment.slice_velocity_range(ranges[2],ranges[3]);
	mid_segments=[];
	slow_segments=[];
	for s in out_segments:
		mid_tmp, slow_tmp = s.slice_velocity_range(ranges[0],ranges[1]);
		mid_segments = mid_segments + mid_tmp;
		slow_segments = slow_segments + slow_tmp;
	js_file.write("var summary_strings = Array();");
	js_file.write("var segments_array = Array();");
	combined_length = 0.0
	combined_points = 0;

	for i in fast_segments + mid_segments + slow_segments:
		combined_points += i.num_points()
	assert combined_points == segment.num_points(), "lost some points (had %d, now %d)" %(segments.num_points(), combined_points)
		
#	sa,sb,sc = segment.slice_largest_constant_velocity_segment(20);
	for s in fast_segments:
		js_file.write("segments_array.push(%s);" % s.javascript_string("#00ff00"));
	for s in mid_segments:
		js_file.write("segments_array.push(%s);" % s.javascript_string("#ffff00"));
	for s in slow_segments:
		js_file.write("segments_array.push(%s);" % s.javascript_string("#ff0000"));

	js_file.write(segment.summary_string("summary_strings.push(\"", "\");\n"));

	# write a corresponding html file 
	html_filename = js_filename.replace("json","html"); 
	output_html = open(html_filename,"w");
	template_html_filename = open("html/data.html.template","r")
	output_html.write(template_html_filename.read().replace("__JSON_FILENAME__",js_filename[5:]))

	
