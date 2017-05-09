#!/bin/R

#........................................................
# Install packages and set connections  -----------------
# 
# Connect to database and ...... 
#
#.........................................................

install.packages('RPostgreSQL')
install.packages('rgrass7')
install.packages('igraph')
install.packages('doMC')

#Setup connection from R to PostGIS
require(RPostgreSQL)

#Set connection parameters
pg_drv<-dbDriver("PostgreSQL")
pg_db <- "nofa_sandbox"
pg_schema <- "Hydrography"
pg_tmp_schema <- "temporary"
pg_user <- "sam.perrin"
pg_password <- "vegemite"
pg_host <- "vm-srv-finstad.vm.ntnu.no"

wrid <- 14695
dem <- "dem_10m_nosefi_float"

#Initialise connection
con<-dbConnect(pg_drv,dbname=pg_db,user=pg_user, password=pg_password,host=pg_host)

#.........................................................
# Import spatial data from DB into PostGIS  -----------
#.........................................................


#Overlay river network by major watersheds / catchments 
#(for limiting connectivity analysis to within catchment connectivity)
#remove  "WHERE vassomr = '123'" in order to run on the entire dataset
#this joins stream IDs with Water Region IDs, only needs to be used if tables need to be updated
dbSendQuery(con,paste("DROP TABLE IF EXISTS ", pg_tmp_schema, ".\"Streams_WR", wrid, "\"", sep=''))
dbSendQuery(con,paste("CREATE TABLE ", pg_tmp_schema, ".\"Streams_WR", wrid, "\" AS SELECT a.*, max(b.gid) AS \"waterRegionID\", array_agg(b.gid) AS \"waterRegionIDs\" FROM
\"Hydrography\".\"Streams_Norway_2017\" AS a,
(SELECT gid, geom FROM \"Hydrography\".\"waterregions_dem_10m_nosefi\" WHERE gid = ", wrid, ") AS b
WHERE ST_Intersects(a.geom,b.geom)
GROUP BY a.ogc_fid, b.gid
ORDER BY a.ogc_fid, b.gid;", sep=''))

#Create spatial index and analyse table in order to speed up subsequent queries 
# (VACUUM ANALYZE needed in order to make PostGIS use the spatial index)
#Only works if you run series of commands directly above this
#Basically the index helps speed things up by letting the connection know if a stream is close to our watershed and worth checking further to see if there is intersection
dbSendQuery(con,paste("CREATE INDEX \"Streams_WR", wrid, "_gist\" ON ", pg_tmp_schema, ".\"Streams_WR", wrid, "\" USING gist (geom);", sep=''))
dbSendQuery(con,paste("ALTER TABLE ", pg_tmp_schema, ".\"Streams_WR", wrid, "\" CLUSTER ON \"Streams_WR", wrid, "_gist\";", sep=''))
dbSendQuery(con,paste("VACUUM FULL ANALYZE ", pg_tmp_schema, ".\"Streams_WR", wrid, "\"", sep=''))

#Same as above but for lakes
#Overlay lakes by major watersheds / catchments 
#(for limiting connectivity analysis to within catchment connectivity)
#remove  "WHERE vassomr = '123'" in order to run on the entire dataset
dbSendQuery(con,paste("DROP TABLE IF EXISTS ", pg_tmp_schema, ".\"Lakes_WR", wrid, "\"", sep=''))
dbSendQuery(con,paste("CREATE TABLE ", pg_tmp_schema, ".\"Lakes_WR", wrid, "\" AS SELECT a.gid AS lake_id, a.geom, a.cat AS lake_gid, a.n, max(b.gid) AS \"waterRegionID\", array_agg(b.gid) AS \"waterRegionIDs\" FROM 
temporary.lakes_nosefi AS a,
(SELECT gid, geom FROM \"Hydrography\".\"waterregions_dem_10m_nosefi\" WHERE gid = ", wrid, ") AS b
WHERE ST_Intersects(a.geom, b.geom)
GROUP BY a.gid, a.geom, a.cat, a.n
ORDER BY a.gid;", sep=''))

#Create spatial index and analyse table in order to speed up subsequent queries (VACUUM ANALYZE needed in order to make PostGIS use the spatial index)
dbSendQuery(con,paste("CREATE INDEX \"Lakes_WR", wrid, "_gist\" ON ", pg_tmp_schema, ".\"Lakes_WR", wrid, "\" USING gist (geom);", sep=''))
dbSendQuery(con,paste("ALTER TABLE ", pg_tmp_schema, ".\"Lakes_WR", wrid, "\" CLUSTER ON \"Lakes_WR", wrid, "_gist\";", sep=''))
dbSendQuery(con,paste("VACUUM FULL ANALYZE ", pg_tmp_schema, ".\"Lakes_WR", wrid, "\"", sep=''))

#Remove interior rings (islands) from lake polygons, because NVEs rivernetwork is crossing islands (which is only a technical artifact)
dbSendQuery(con, paste("DROP TABLE IF EXISTS ", pg_tmp_schema, ".\"Lakes_WR", wrid, "_erings\"", sep=''))
dbSendQuery(con, paste("CREATE TABLE ", pg_tmp_schema, ".\"Lakes_WR", wrid, "_erings\" AS SELECT lake_id, lake_gid, geom, n FROM
 (SELECT lake_id, lake_gid, n, (ST_Dump(ST_Union(ST_BuildArea(ST_ExteriorRing(geom))))).geom AS geom FROM 
 (SELECT lake_id, lake_gid, n, (ST_Dump(geom)).geom AS geom FROM ", pg_tmp_schema, ".\"Lakes_WR", wrid, "\") AS a GROUP BY lake_id, lake_gid, n) AS b;", sep=''))

#Create spatial index and analyse table in order to speed up subsequent queries (VACUUM ANALYZE needed in order to make PostGIS use the spatial index)
dbSendQuery(con, paste("CREATE INDEX \"Lakes_WR", wrid, "_erings_gist\" ON ", pg_tmp_schema, ".\"Lakes_WR", wrid, "_erings\" USING gist (geom);", sep=''))
dbSendQuery(con, paste("ALTER TABLE ", pg_tmp_schema, ".\"Lakes_WR", wrid, "_erings\" CLUSTER ON \"Lakes_WR", wrid, "_erings_gist\";", sep=''))
dbSendQuery(con, paste("VACUUM FULL ANALYZE ", pg_tmp_schema, ".\"Lakes_WR", wrid, "_erings\"", sep=''))

# #Union lake polygons by individual ID in order to use the "real" lake polygon for displaying purposes
# #Each lake should only occure in 1 major watershed and should have only 1 "vatnlnr" (because of the above query) possible numbering errors (e.g. lakes consisting of 2 polygons with different "vatnlnr"s are documented in column "vatnlnrs", an array containg all vatnlnrs of the unioned lake polygon)
# dbSendQuery(con, paste("DROP TABLE IF EXISTS kladd.n50_2013_innsjoe_uniq_vassdragsomraader;", sep=''))
# dbSendQuery(con, paste("CREATE TABLE kladd.n50_2013_innsjoe_uniq_vassdragsomraader AS SELECT b.lake_id, max(vatnlnr) AS vatnlnr, uniq(sort(array_agg(vatnlnr))) AS vatnlnrs, array_length(uniq(sort(array_agg(vatnlnr))), 1) AS vatnlnr_c, vassomr, ST_Union(a.geom) AS geom FROM
# kladd.n50_2013_innsjoe_vassdragsomraader AS a,
# kladd.n50_2013_innsjoe_vassdragsomraader_erings AS b
# WHERE ST_Intersects(a.geom,b.geom)
# GROUP BY b.lake_id, vassomr;", sep=''))

# #Create spatial index and analyse table in order to speed up subsequent queries (VACUUM ANALYZE needed in order to make PostGIS use the spatial index)
# dbSendQuery(con, paste("CREATE INDEX n50_2013_innsjoe_uniq_vassdragsomraader_gist ON kladd.n50_2013_innsjoe_uniq_vassdragsomraader USING gist (geom);", sep=''))
# dbSendQuery(con, paste("ALTER TABLE kladd.n50_2013_innsjoe_uniq_vassdragsomraader CLUSTER ON n50_2013_innsjoe_uniq_vassdragsomraader_gist;", sep=''))
# dbSendQuery(con, paste("VACUUM FULL ANALYZE kladd.n50_2013_innsjoe_uniq_vassdragsomraader", sep=''))

#Disconnect from PostGIS
dbDisconnect(con)

#.........................................................
# Build network in GRASS GIS  ------
#.........................................................


##Connect to GRASS
library("rgrass7")
#create a new mapset (use only if necessary)
#try(system("grass72 -text -c -e /media/harddisk/grassdata/ETRS_33N/p_INVAFISH_nidelva"))
initGRASS(gisBase ='/usr/local/grass-7.2.1svn/', location = 'ETRS_33N', mapset = 'u_sam.perrin',
          gisDbase = '/media/harddisk/grassdata', override = TRUE)

#Import relevant tables from PostGIS to GRASS (v.external as a possible and faster alternative?)
# "overwrite" gets rid of tables if they are already there so they can be updated, "verbose" brings in extra detail
#Username and password!!!
execGRASS("v.in.ogr", flags=c("o", "overwrite", "verbose"), 
          input=paste("PG:dbname=", pg_db, "active_schema=", pg_tmp_schema, sep=''), 
		  layer=paste("\"Streams_WR", wrid, "\"", sep=''), 
		  output=paste("Streams_WR", wrid, sep=''), snap=0.001)
execGRASS("v.in.ogr", flags=c("o", "overwrite", "verbose"),
          input=paste("PG:dbname=", pg_db, "active_schema=", pg_tmp_schema, sep=''),
		  layer=paste("\"Lakes_WR", wrid, "_erings\"", sep=''),
		  output=paste("Lakes_WR", wrid, "_erings", sep=''), snap=0.001)
execGRASS("v.in.ogr", flags=c("o", "overwrite", "verbose"),
          input=paste("PG:dbname=", pg_db, "active_schema=", pg_tmp_schema, sep=''),
		  layer=paste("\"Lakes_WR", wrid, "\"", sep=''),
		  output=paste("Lakes_WR", wrid, sep=''), snap=0.001)

#Merge adjacent lines as far as possible (in order to have only single line strings between lakes / fork points (lake - lake, lake - fork point, fork point - fork point))
#If you have a stream that has been broken into two sections by a rogue node somewhere, merges them into 1
execGRASS("v.build.polylines", flags=c("overwrite", "verbose"),
          input=paste("Streams_WR", wrid, sep=''),
		  output=paste("Streams_WR", wrid, "_merge", sep=''), cats="first")
#Remove lines within lakes from river network
execGRASS("v.overlay", flags=c("overwrite", "verbose"),
          ainput=paste("Streams_WR", wrid, "_merge", sep=''), atype="line",
		  binput=paste("Lakes_WR", wrid, "_erings", sep=''),
		  output=paste("Streams_WR", wrid, "_notLake", sep=''),
		  operator="not", olayer="1,1,0")
#Snap nodes in rivernetwork (in order to reduce number of topological errors)
#If for some reason you have a stream that is broken up, this snaps them together
execGRASS("v.edit", map=paste("Streams_WR", wrid, "_notLake", sep=''),
          type="line", tool="snap", threshold=c(-1,1,0),
          where="cat>0", snap="node")

#This MIGHT snap rivers to lakes
execGRASS("v.edit", map=paste("Streams_WR", wrid, "_notLake", sep=''),
          type="line", tool="snap", threshold=c(-1,1,0),
          where="cat>0", snap="vertex")

#Create a network dataset from river line strings by adding nodes at start, end, and fork points
execGRASS("v.net", flags=c("c","overwrite"),
          input=paste("Streams_WR", wrid, "_notLake", sep=''),
		  output=paste("Streams_WR", wrid, "_network", sep=''),
		  operation="nodes")


#Assign ID of closest lake to nodes in network
#create a database table for vertices
execGRASS("v.db.addtable", flags=c("verbose"),
          map=paste("Streams_WR", wrid, "_network", sep=''), layer=2,
		  table=paste("Streams_WR", wrid, "_network_vertices", sep=''),
		  columns="lake_id integer")
#Load ID of lakes within 1m distance to vertices into database of vertices (column = lake_id)
execGRASS("v.distance", flags=c("overwrite","verbose"),
          from=paste("Streams_WR", wrid, "_network", sep=''), from_layer="2", column="lake_id",
		  to=paste("Lakes_WR", wrid, "_erings", sep=''), to_column="lake_id",
		  upload="to_attr", dmax=1)
##Mark lake vertices
#execGRASS("v.db.update", map=paste("Streams_WR", wrid, "_network", sep=''), layer="2", column="is_lake", value=1, where="lake_id IS NULL",redirect=TRUE,legacyExec=TRUE)

#Fill lake_id column for vertices with no lake within 1m distance with individual values from "cat" column
max_lake_id <- max(as.integer(as.character(execGRASS("v.db.select", flags=c("c", "quiet"),map=paste("Lakes_WR", wrid, "_erings", sep=''), layer="1", columns="lake_id", redirect=TRUE, legacyExec=TRUE))), na.omit=TRUE)


execGRASS("v.db.update",
          map=paste("Streams_WR", wrid, "_network", sep=''), layer="2", column="lake_id",
		  query_column=paste("cat+", max_lake_id, sep=""),
		  where="lake_id IS NULL",redirect=TRUE,legacyExec=TRUE)
		  
#.........................................................
# Calculate stream/lake data based on spatial data -----
#.........................................................


#Calculate lake size
execGRASS("v.db.addcolumn",
          map=paste("Lakes_WR", wrid, sep=''),
		  columns="area_ha double precision")
execGRASS("v.to.db", flags="quiet",
          map=paste("Lakes_WR", wrid, sep=''),
		  option="area", columns="area_ha", units="hectares")
lake_feat <- data.frame(matrix(unlist(strsplit(execGRASS("v.db.select", flags=c("c", "quiet"),
                                                         map=paste("Lakes_WR", wrid, sep=''),
														 layer="1", columns="lake_id,area_ha",
														 redirect=TRUE, separator=" ", legacyExec=TRUE
														 ), " ")),
														 nrow=vDataCount(paste("Lakes_WR", wrid, sep=''), layer="1"),
														 byrow=T))
names(lake_feat) <- c("lake_id","area_ha")

##Charaterize river line strings
#Convert line strings from vector to raster using ID ("cat")
execGRASS("g.region", flags="p",
          vector=c(paste("Lakes_WR", wrid, "_erings", sep=''),
		           paste("Streams_WR", wrid, "_network", sep='')),
		  n="n+100", s="s-100", e="e+100", w="w-100",
		  align=dem)
execGRASS("v.to.rast", flags=c("overwrite", "verbose"),
          input=paste("Streams_WR", wrid, "_network", sep=''), type="line",
		  output=paste("Streams_WR", wrid, "_network_cat", sep=''), use="cat")

#Calculate slope in river network
#First we look at the pixels directly lef,right,above,below central pixel, see if the river is present there, 
#take the altitude change, divide by 10(metres)
execGRASS("r.mapcalc", flags="overwrite",
          expression=paste("Streams_WR", wrid, "_network_local_slope_direct=if(isnull(Streams_WR", 
		                   wrid, "_network_cat),null(),if(nmin(if(isnull(Streams_WR",
						   wrid, "_network_cat[1,0]),9999,", dem, "[1,0]),if(isnull(Streams_WR",
						   wrid, "_network_cat[-1,0]),9999,", dem, "[-1,0]),if(isnull(Streams_WR",
						   wrid, "_network_cat[0,1]),9999,", dem, "[0,1]),if(isnull(Streams_WR",
						   wrid, "_network_cat[0,-1]),9999,", dem, "[0,-1]))>", dem, ",-9999,(", dem, "-nmin(if(isnull(Streams_WR",
						   wrid, "_network_cat[1,0]),9999,", dem, "[1,0]),if(isnull(Streams_WR",
						   wrid, "_network_cat[-1,0]),9999,", dem, "[-1,0]),if(isnull(Streams_WR",
						   wrid, "_network_cat[0,1]),9999,", dem, "[0,1]),if(isnull(Streams_WR",
						   wrid, "_network_cat[0,-1]),9999,", dem, "[0,-1])))/10.0))", sep=''))
#then take pixels diagonal from central pixel, divide by sqrt(200)
execGRASS("r.mapcalc", flags="overwrite",
          expression=paste("Streams_WR", wrid, "_network_local_slope_diagonal=if(isnull(Streams_WR", 
		                   wrid, "_network_cat),null(),if(nmin(if(isnull(Streams_WR", 
		                   wrid, "_network_cat[1,1]),9999,", dem, "[1,1]),if(isnull(Streams_WR", 
		                   wrid, "_network_cat[-1,1]),9999,", dem, "[-1,1]),if(isnull(Streams_WR", 
		                   wrid, "_network_cat[1,-1]),9999,", dem, "[1,-1]),if(isnull(Streams_WR", 
		                   wrid, "_network_cat[-1,-1]),9999,", dem, "[-1,-1]))>", dem, ",-9999,(", dem, "-nmin(if(isnull(Streams_WR", 
		                   wrid, "_network_cat[1,1]),9999,", dem, "[1,1]),if(isnull(Streams_WR", 
		                   wrid, "_network_cat[-1,1]),9999,", dem, "[-1,1]),if(isnull(Streams_WR", 
		                   wrid, "_network_cat[1,-1]),9999,", dem, "[1,-1]),if(isnull(Streams_WR", 
		                   wrid, "_network_cat[-1,-1]),9999,", dem, "[-1,-1])))/sqrt(200.0)))", sep=''))
#now take the maximum of all of them
execGRASS("r.mapcalc", flags="overwrite",
          expression=paste("Streams_WR", wrid, "_network_local_slope=max(if(Streams_WR",
		                   wrid, "_network_local_slope_direct==-9999,0,Streams_WR",
						   wrid, "_network_local_slope_direct),if(Streams_WR",
						   wrid, "_network_local_slope_diagonal==-9999,0,Streams_WR",
						   wrid, "_network_local_slope_diagonal))", sep=''))
execGRASS("g.remove", type="raster", pattern=paste("Streams_WR", wrid, "_network_local_slope_di*", sep=''), flags="f")

#Calculate slope statistics for river network
#execGRASS("v.rast.stats", flags="c", map="river_lake_network_70m_cat", layer=3, raster="river_lake_network_local_slope", column_prefix="slope")
#execGRASS("v.rast.stats", flags="c", map=paste("Streams_WR", wrid, "_network", sep=''), layer=1, raster="river_lake_network_local_slope", column_prefix="slope")
#execGRASS("v.rast.stats", flags="c", map=paste("Streams_WR", wrid, "_network", sep=''), raster="dem_10m_nosefi_float", column_prefix="altitude")
#Calculate univariate statistics on terrain model for every line string ("cat" as key column)
#First converts line into series of rasters, then gives univar stats for that series
dem_stats_pre <- execGRASS("r.univar", flags=c("t","e","quiet"),
                           map=dem, zones=paste("Streams_WR", wrid, "_network_cat", sep=''),
						   separator=',', percentile=c(10,90), redirect=TRUE, legacyExec=TRUE)
dem_stats <- data.frame(matrix(unlist(strsplit(dem_stats_pre[2:length(dem_stats_pre)], ",")), nrow=length(dem_stats_pre)-1, byrow=T))
names(dem_stats) <- paste("dem_", strsplit(dem_stats_pre[1], ",")[[1]], sep="")
names(dem_stats)[1] <- "line_category"

#Calculate univariate statistics on slope for every line string ("cat" as key column)
slope_stats_pre <- execGRASS("r.univar", flags=c("t","e","quiet"),
                             map=paste("Streams_WR", wrid, "_network_local_slope", sep=''),
                             zones=paste("Streams_WR", wrid, "_network_cat", sep=''),
							 separator=',', percentile=c(10,90), redirect=TRUE, legacyExec=TRUE)

slope_stats <- data.frame(matrix(unlist(strsplit(slope_stats_pre[2:length(slope_stats_pre)], ",")), nrow=length(slope_stats_pre)-1, byrow=T))
names(slope_stats) <- paste("slope_", strsplit(slope_stats_pre[1], ",")[[1]], sep="")
names(slope_stats)[1] <- "line_category"

#Calculate lenght in meter for every line string ("cat" as key column)
execGRASS("v.db.addcolumn",
          map=paste("Streams_WR", wrid, "_network", sep=''),
		  columns="length_m double precision")
execGRASS("v.to.db", flags="quiet",
          map=paste("Streams_WR", wrid, "_network", sep=''), type="line",
		  option="length", columns="length_m", units="meters")

line_length <- data.frame(matrix(unlist(strsplit(execGRASS("v.db.select", flags=c("c", "quiet"),
                                                           map=paste("Streams_WR", wrid, "_network", sep=''), layer="1",
														   columns="cat,length_m", redirect=TRUE,
														   separator=" ", legacyExec=TRUE), " ")),
										nrow=vDataCount(paste("Streams_WR", wrid, "_network", sep=''), layer="1"), byrow=T))
names(line_length) <- c("line_category", "length_m")


#.........................................................
## Complete network analysis using igraph package -------
#.........................................................

##Load network data to R
#Vertices
v <- data.frame(matrix(unlist(strsplit(execGRASS("v.db.select", flags=c("c", "quiet"),
                                                 map=paste("Streams_WR", wrid, "_network", sep=''), layer="2",
												 redirect=TRUE, separator=" ", legacyExec=TRUE), " ")),
							nrow=vDataCount(paste("Streams_WR", wrid, "_network", sep=''), layer="2"), byrow=T))
names(v) <- c("point_category", "lake_id")
v_uniq <- data.frame(lake_id=unique(as.character(v$lake_id)))
v_uniq <- merge(v_uniq, lake_feat, all.x=TRUE)
v_uniq$is_lake=ifelse(as.numeric(as.character(v_uniq$area_ha))>0, 1, 0)

v_uniq$lake_area_ha[is.na(v_uniq$area_ha)] <- 0
v_uniq$is_lake[is.na(v_uniq$is_lake)] <- 0


#Edges
e <- data.frame(matrix(unlist(strsplit(execGRASS("v.net",
                                                 input=paste("Streams_WR", wrid, "_network", sep=''),
												 points=paste("Streams_WR", wrid, "_network", sep=''),
												 operation="report",
												 arc_layer="1", node_layer="2",
												 flags="quiet", redirect=TRUE, legacyExec=TRUE), " ")),
							nrow=vDataCount(paste("Streams_WR", wrid, "_network", sep=''), layer="1"), byrow=T))
names(e) <- c("line_category","start_point_category","end_point_category")
#e <- merge(merge(e,data.frame(start_point_category=v$point_category, from_lake=v$lake_id)),data.frame(end_point_category=v$point_category, to_lake=v$lake_id))
e <- merge(merge(e, data.frame(end_point_category=v$point_category, to_lake_id=v$lake_id)),
                 data.frame(start_point_category=v$point_category, from_lake_id=v$lake_id))
#Remove loops (edges between identical vertices)
e <- e[grep(TRUE, e$from_lake_id != e$to_lake_id),]

#e_attributes <- data.frame(matrix(unlist(strsplit(execGRASS("v.db.select", flags="c", map=paste("Streams_WR", wrid, "_network", sep=''), columns=c("cat,altitude_average,slope_average,length_m"), layer="1", redirect=TRUE, separator=" ", legacyExec=TRUE), " ")), nrow=vDataCount(paste("Streams_WR", wrid, "_network", sep=''), layer="1"), byrow=T))
#names(e_attributes) <- c("line_category","altitude_average","slope_average","length_m")

e <- merge(merge(merge(e,line_length),dem_stats),slope_stats)
#e <- merge(e,e_attributes)

# Classify nodes as 1 = Sources (without incomming edges), 2 = Connections (one incoming and one outgoing edge), 3 = Fork (two or more incoming edges), 4 = Outlets (without outgoing edges)
v_uniq$typ <- 0 #Unclassified
v_uniq$typ[(as.character(v_uniq$lake_id) %in% unique(sort(as.character(e$to_lake_id)))) == FALSE] <- 1 # Sources (without incomming edges)
v_uniq$typ[(as.character(v_uniq$lake_id) %in% unique(sort(as.character(e$from_lake_id)))) == FALSE] <- 4 # Outlets (without outgoing edges)
v_uniq$typ[(as.character(v_uniq$lake_id) %in% unique(sort(as.character(e$from_lake_id)))) == FALSE && 
           (as.character(v_uniq$lake_id) %in% unique(sort(as.character(e$to_lake_id)))) == FALSE] <- 5 # Isolated nodes (without incoming or outgoing edges)
#v_uniq$typ[(as.character(v_uniq$lake_id) %in% unique(sort(as.character(e$from_lake_id)))) == FALSE && 
#           (as.character(v_uniq$lake_id) %in% unique(sort(as.character(e$to_lake_id)))) == FALSE] <- 5 # Isolated nodes (without incoming or outgoing edges)
v_uniq$typ <- ifelse(v_uniq$typ==0,ifelse(table(as.character(e$to_lake_id))>1,3,2),v_uniq$typ)

library(igraph)
library(doMC)
registerDoMC()

##Build graph
#Create empty, directed graph
g <- graph.empty()
#Add vertices
g <- add.vertices(g, nrow(v_uniq),
                  lake_id=as.integer(as.character(v_uniq$lake_id)),
				  is_lake=as.integer(as.character(v_uniq$is_lake)),
				  lake_area_ha=as.double(as.character(v_uniq$area_ha)),
				  typ=as.integer(as.character(v_uniq$typ)))
#Define edges
names <- V(g)$lake_id
ids <- 1:length(names)
names(ids) <- names
from <- as.character(e$from_lake_id)
to <- as.character(e$to_lake_id)
edges <- matrix(c(ids[from], ids[to]), nc=2)

#Add edges in downstream direction to graph
g_down <- add.edges(g, t(edges), cat=e$line_category, from=e$from_lake_id, to=e$to_lake_id,
                    length_m=as.double(as.character(e$length_m)),
					altitude_average=as.double(as.character(e$dem_mean)),
					altitude_min=as.double(as.character(e$dem_min)),
					altitude_max=as.double(as.character(e$dem_max)),
					altitude_stddev=as.double(as.character(e$dem_stddev)),
					slope_min=as.double(as.character(e$slope_min)),
					slope_perc_10=as.double(as.character(e$slope_perc_10)),
					slope_first_quart=as.double(as.character(e$slope_first_quart)),
					slope_mean=as.double(as.character(e$slope_mean)),
					slope_third_quart=as.double(as.character(e$slope_third_quart)),
					slope_perc_90=as.double(as.character(e$slope_perc_90)),
					slope_max=as.double(as.character(e$slope_max)),
					slope_stddev=as.double(as.character(e$slope_stddev)),
					slope_variance=as.double(as.character(e$slope_variance)),
					direction="downstreams")

#Add edges in upstream direction to graph
edges_u <- matrix(c(ids[to], ids[from]), nc=2)
g_up <- add.edges(g, t(edges_u), cat=e$line_category, from=e$from_lake_id, to=e$to_lake_id, 
                  length_m=as.double(as.character(e$length_m)),
				  altitude_average=as.double(as.character(e$dem_mean)),
				  altitude_min=as.double(as.character(e$dem_min)),
				  altitude_max=as.double(as.character(e$dem_max)),
				  altitude_stddev=as.double(as.character(e$dem_stddev)),
				  slope_min=as.double(as.character(e$slope_min)),
				  slope_perc_10=as.double(as.character(e$slope_perc_10)),
				  slope_first_quart=as.double(as.character(e$slope_first_quart)),
				  slope_mean=as.double(as.character(e$slope_mean)),
				  slope_third_quart=as.double(as.character(e$slope_third_quart)),
				  slope_perc_90=as.double(as.character(e$slope_perc_90)),
				  slope_max=as.double(as.character(e$slope_max)),
				  slope_stddev=as.double(as.character(e$slope_stddev)),
				  slope_variance=as.double(as.character(e$slope_variance)),
				  direction="upstreams")

#Add edges in bi-directional graph
g <- add.edges(g_down, t(edges_u), cat=e$line_category, from=e$from_lake_id, to=e$to_lake_id, 
               length_m=e$length_m,
			   length_m=as.double(as.character(e$length_m)),
			   altitude_average=as.double(as.character(e$dem_mean)),
			   altitude_min=as.double(as.character(e$dem_min)),
			   altitude_max=as.double(as.character(e$dem_max)),
			   altitude_stddev=as.double(as.character(e$dem_stddev)),
			   slope_min=as.double(as.character(e$slope_min)),
			   slope_perc_10=as.double(as.character(e$slope_perc_10)),
			   slope_first_quart=as.double(as.character(e$slope_first_quart)),
			   slope_mean=as.double(as.character(e$slope_mean)),
			   slope_third_quart=as.double(as.character(e$slope_third_quart)),
			   slope_perc_90=as.double(as.character(e$slope_perc_90)),
			   slope_max=as.double(as.character(e$slope_max)),
			   slope_stddev=as.double(as.character(e$slope_stddev)),
			   slope_variance=as.double(as.character(e$slope_variance)),
			   direction="bi-directional")

#V(g)$typ <- 0 #Unclassified
#V(g)$typ[(as.character(V(g)$lake_id) %in% as.character(E(g)$to)) == FALSE] <- 1 # Sources (without incomming edges)
#V(g)$typ[(as.character(V(g)$lake_id) %in% as.character(E(g)$from)) == FALSE] <- 2 # Outlets (without outgoing edges)
#V(g)$typ <- ifelse(V(g)$typ==0,ifelse(table(E(g)$to)>2,3,4),V(g)$typ) # 3 = Connections (one incoming and one outgoing edge), 4 = Fork (two or more incoming edges)

#######################################################################
#Reduce number of possible connections:
# - Use only upstreams connections (drop all downstreams combinations)
# - Use only direct / nearest neighbor connections
# - Use only direct / nearest neighbor connections upstreams
#
#Some simple connectivity measures:
# - Clustersize (in number lakes)
# - Clustersize (in lake area)
# - length of river network downstreams
# - number of lakes downstreams
# - area of lakes downstreams
# - list of lakes downstreams
# - length of river network upstreams
# - number of lakes upstreams
# - area of lakes upstreams
# - list of lakes upstreams:
#   1 2,9,10,11
#   2 9,10,11
#   3 4,5,16
#   4 5,16
# - list of first lakes upstreams:
#   1 2
#   2 9
#   3 4
#   4 5, 16
#
# - connectivity to first lakes upstreams, characterised by:
#   - total length
#   - min slope (average over all river sections)
#   - 10% percentile of slope (average over all river sections)
#   - first quartile of slope (average over all river sections)
#   - mean slope (average over all river sections)
#   - third quartile of slope (average over all river sections)
#   - 90% percentile of slope (average over all river sections)
#   - max slope (average over all river sections)
##   - length of river sections by slope classes (slope < X %, slope >= X % AND slope < Y %, ...)
##   - number of river sections with slope > X % (a bit more complicated and probably less reliable)
##   - max altitude
##   - contains dam(s), whire(s) (boolean)
##   - contained dams, whires (list of contained dams and whires)
## Stepwise analysis
## direct neighbors and downsreams sequence
## habitat models (per species)
## likelyhood to get from Lake A to B
##   - use lakes with invasion
## lieklyhood to survive and reproduce in Lake X or Y


#.........................................................
# Create network vaiables using R -------


cores <- 5

require(postGIStools)
lakes_withData <- 'SELECT * FROM temporary."lakeID_withData" '
lakes_withDataT <- get_postgis_query(con,lakes_withData)
lakes <- lakes_withDataT$lake_id
lakes <- lakes[!is.na(lakes)]


lakes <- V(g_down)$lake_id[V(g_down)$is_lake==1]

str(lakes)


cl <- clusters(g_down)
V(g_down)$cluster <- cl$membership
V(g_down)$cluster_lake_area_ha <- unlist(mclapply(1:length(V(g_down)$cluster), function(c) sum(V(g_down)$lake_area_ha[V(g_down)$cluster==V(g_down)$cluster[c]],na.rm=TRUE), mc.cores=cores))
V(g_down)$cluster_lake_n <- unlist(mclapply(1:length(V(g_down)$cluster), function(c) sum(V(g_down)$is_lake[V(g_down)$cluster==V(g_down)$cluster[c]]), mc.cores=cores))

downstream_paths <- lapply(1:length(lakes), function(l) get.shortest.paths(g_down, 
                from=V(g_down)[V(g_down)$lake_id==lakes[l]], 
                to=V(g_down)[grep(TRUE, V(g_down)$typ==4 &
                V(g_down)$cluster==V(g_down)$cluster[V(g_down)$lake_id==lakes[l]])],
                mode=c("out"), output="vpath")[[1]][[1]])




# calculate for every lake the path to the outlet (as vpath)
downstream_paths <- mclapply(1:length(lakes), function(l) get.shortest.paths(g_down, from=V(g_down)[V(g_down)$lake_id==lakes[l]], to=V(g_down)[grep(TRUE, V(g_down)$typ==4 & V(g_down)$cluster==V(g_down)$cluster[V(g_down)$lake_id==lakes[l]])], mode=c("out"), output="vpath")[[1]][[1]], mc.cores = cores)
# remove all non-lake vertices from vpath (from-lake is allways included in the resulting path and has to be removed)
downstream_lakes <- mclapply(1:length(lakes), function(l) as.integer(grep(lakes[l], V(g_down)[unlist(downstream_paths[l])[V(g_down)[unlist(downstream_paths[l])]$is_lake==1]]$lake_id, value=TRUE, invert=TRUE)), mc.cores = cores)
# concatenate lake_ids to string
downstream_lakes_str <- mclapply(1:length(lakes), function(l) toString(downstream_lakes[[l]]), mc.cores = cores)
# calculate number of lakes downstreams
downstream_lakes_n <- mclapply(1:length(lakes), function(l) ifelse(length(downstream_lakes[[l]])==0,0,length(downstream_lakes[[l]])-1), mc.cores = cores)
# calculate area of lakes downstreams
downstream_lakes_area_ha <- mclapply(1:length(lakes), function(l) sum(V(g_down)[unlist(downstream_paths[l])[V(g_down)[unlist(downstream_paths[l])]$is_lake==1]]$lake_area_ha), mc.cores = cores)
# get first lake downstreams
first_downstream_lake <- mclapply(1:length(lakes), function(l) downstream_lakes[[l]][1], mc.cores = cores)


downstream_epaths <- mclapply(1:length(lakes), function(l) get.shortest.paths(g_down, from=V(g_down)[V(g_down)$lake_id==lakes[l]], to=V(g_down)[grep(TRUE, V(g_down)$typ==4 & V(g_down)$cluster==V(g_down)$cluster[V(g_down)$lake_id==lakes[l]])], mode=c("out"), output="epath")$epath[[1]], mc.cores = cores)
downstream_stream_length_km <- mclapply(1:length(lakes), function(l) sum(E(g_down)[downstream_epaths[[l]]]$length_m)/1000.0, mc.cores = cores)


cl <- clusters(g_up)
V(g_up)$cluster <- cl$membership

# calculate for every lake the path to the sources (as vpath)
upstream_paths <- mclapply(1:length(lakes), function(l) unique(sort(unlist(get.shortest.paths(g_up, from=V(g_up)[V(g_up)$lake_id==lakes[l]], to=V(g_up)[grep(TRUE, V(g_up)$typ==1 & V(g_up)$cluster==V(g_up)$cluster[V(g_up)$lake_id==lakes[l]])], mode=c("out"), output="vpath")$vpath))), mc.cores = cores)
# remove all non-lake vertices from vpath
upstream_lakes <- mclapply(1:length(lakes), function(l) as.integer(grep(lakes[l], V(g_up)[upstream_paths[[l]][V(g_up)[upstream_paths[[l]]]$is_lake==1]]$lake_id, value=TRUE, invert=TRUE)), mc.cores = cores)
# concatenate lake_ids to string
upstream_lakes_str <- mclapply(1:length(lakes), function(l) toString(upstream_lakes[[l]]), mc.cores = cores)
# calculate number of lakes upstreams
upstream_lakes_n <- mclapply(1:length(lakes), function(l) length(upstream_lakes[[l]]), mc.cores = cores)
# calculate area of lakes upstreams
upstream_lakes_area_ha <- mclapply(1:length(lakes), function(l) sum(V(g_up)$lake_area_ha[grep(TRUE, V(g_up)$lake_id %in% upstream_lakes[[l]])]), mc.cores = cores)

upstream_epaths <- mclapply(1:length(lakes), function(l) unique(sort(unlist(get.shortest.paths(g_up, from=V(g_up)[V(g_up)$lake_id==lakes[l]], to=V(g_up)[grep(TRUE, V(g_up)$typ==1 & V(g_up)$cluster==V(g_up)$cluster[V(g_up)$lake_id==lakes[l]])], mode=c("out"), output="epath")$epath))), mc.cores = cores)

upstream_stream_length_km <- mclapply(1:length(lakes), function(l) sum(E(g_up)[upstream_epaths[[l]]]$length_m)/1000.0, mc.cores = cores)

# get first lake upstreams

first_upstream_lake_paths <- mclapply(1:length(lakes), function(l) get.shortest.paths(g_up,
                                                                                      from=V(g_up)[grep(TRUE, V(g_up)$lake_id==lakes[l])],
                                                                                      to=V(g_up)[grep(TRUE, V(g_up)$lake_id %in% upstream_lakes[[l]])],
                                                                                      mode=c("out"), output="vpath")$vpath, mc.cores = cores)

first_upstream_lakes <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lake_paths[[l]])>0) {unlist(lapply(1:length(first_upstream_lake_paths[[l]]), function(x) if(length(V(g_up)[unlist(first_upstream_lake_paths[[l]][[x]])[V(g_up)[unlist(first_upstream_lake_paths[[l]][[x]])]$is_lake==1]]$lake_id)==2) {V(g_up)[unlist(first_upstream_lake_paths[[l]][[x]])[V(g_up)[unlist(first_upstream_lake_paths[[l]][[x]])]$is_lake==1]]$lake_id[2]}))}, mc.cores = cores)

first_upstream_lakes_str <- mclapply(1:length(lakes), function(l) toString(first_upstream_lakes[[l]]), mc.cores = cores)


first_upstream_lakes_epath <- mclapply(1:length(lakes), function(l) get.shortest.paths(g_up, from=V(g_up)[V(g_up)$lake_id==lakes[l]], to=V(g_up)[grep(TRUE, V(g_up)$lake_id %in% first_upstream_lakes[[l]])], mode=c("out"), output="epath")$epath, mc.cores = cores)

first_upstream_lakes_slope_min <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_min))}, mc.cores = cores)

first_upstream_lakes_slope_min_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_min[[l]])), mc.cores = cores)

first_upstream_lakes_slope_perc_10 <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_perc_10))}, mc.cores = cores)

first_upstream_lakes_slope_perc_10_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_perc_10[[l]])), mc.cores = cores)

first_upstream_lakes_slope_first_quart <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_first_quart))}, mc.cores = cores)

first_upstream_lakes_slope_first_quart_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_first_quart[[l]])), mc.cores = cores)

first_upstream_lakes_slope_mean <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_mean))}, mc.cores = cores)

first_upstream_lakes_slope_mean_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_mean[[l]])), mc.cores = cores)

first_upstream_lakes_slope_third_quart <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_third_quart))}, mc.cores = cores)

first_upstream_lakes_slope_third_quart_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_third_quart[[l]])), mc.cores = cores)

first_upstream_lakes_slope_perc_90 <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_perc_90))}, mc.cores = cores)

first_upstream_lakes_slope_perc_90_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_perc_90[[l]])), mc.cores = cores)

first_upstream_lakes_slope_max <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_max))}, mc.cores = cores)

first_upstream_lakes_slope_max_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_max[[l]])), mc.cores = cores)

first_upstream_lakes_slope_stddev <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_stddev))}, mc.cores = cores)

first_upstream_lakes_slope_stddev_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_stddev[[l]])), mc.cores = cores)

first_upstream_lakes_slope_variance <- mclapply(1:length(lakes), function(l) if(length(first_upstream_lakes_epath[[l]])>0) {lapply(1:length(first_upstream_lakes_epath[[l]]), function(e) mean(E(g_up)[first_upstream_lakes_epath[[l]][[e]]]$slope_variance))}, mc.cores = cores)

first_upstream_lakes_slope_variance_str <- mclapply(1:length(lakes), function(l) toString(unlist(first_upstream_lakes_slope_variance[[l]])), mc.cores = cores)

# Write data back to PostGIS  -----


result <- data.frame(
lake_id=lakes,
cluster=V(g_down)$cluster[V(g_down)$lake_id %in% lakes],
cluster_lake_area_ha=V(g_down)$cluster_lake_area_ha[V(g_down)$lake_id %in% lakes],
cluster_lake_n=V(g_down)$cluster_lake_n[V(g_down)$lake_id %in% lakes],
downstream_lakes=unlist(downstream_lakes_str),
downstream_lakes_n=unlist(downstream_lakes_n),
downstream_lakes_area_ha=unlist(downstream_lakes_area_ha),
first_downstream_lake=unlist(first_downstream_lake),
downstream_stream_length_km=unlist(downstream_stream_length_km),
upstream_lakes=unlist(upstream_lakes_str),
upstream_lakes_n=unlist(upstream_lakes_n),
upstream_lakes_area_ha=unlist(upstream_lakes_area_ha),
upstream_stream_length_km=unlist(upstream_stream_length_km),
first_upstream_lakes=unlist(first_upstream_lakes_str),
first_upstream_lakes_slope_min=unlist(first_upstream_lakes_slope_min_str),
first_upstream_lakes_slope_perc_10=unlist(first_upstream_lakes_slope_perc_10_str),
first_upstream_lakes_slope_first_quart=unlist(first_upstream_lakes_slope_first_quart_str),
first_upstream_lakes_slope_mean=unlist(first_upstream_lakes_slope_mean_str),
first_upstream_lakes_slope_third_quart=unlist(first_upstream_lakes_slope_third_quart_str),
first_upstream_lakes_slope_perc_90=unlist(first_upstream_lakes_slope_perc_90_str),
first_upstream_lakes_slope_max=unlist(first_upstream_lakes_slope_max_str),
first_upstream_lakes_slope_stddev=unlist(first_upstream_lakes_slope_stddev_str),
first_upstream_lakes_slope_variance=unlist(first_upstream_lakes_slope_variance_str)
)

con<-dbConnect(pg_drv,dbname=pg_db,user=pg_user, password=pg_password,host=pg_host)

result_table <- paste("connectivity_", wrid, sep='')
dbGetQuery(con, paste("SET search_path TO ", pg_tmp_schema, sep=""))
if(dbExistsTable(con, result_table)) {dbRemoveTable(con, result_table)}
dbWriteTable(con, result_table, result)

dbSendQuery(con,paste("CREATE INDEX ", result_table, "_idx ON ", pg_tmp_schema, ".", result_table, " USING btree(lake_id);", sep=''))
dbSendQuery(con,paste("ALTER TABLE ", pg_tmp_schema, ".", result_table, " CLUSTER ON ", result_table, "_idx;", sep=''))
dbSendQuery(con,paste("ALTER TABLE ", pg_tmp_schema, ".", result_table, " ADD CONSTRAINT ", result_table, "_pkey PRIMARY KEY (lake_id);", sep=""))
dbSendQuery(con,paste("VACUUM FULL ANALYZE ", pg_tmp_schema, ".", result_table, ";", sep=''))



execGRASS('v.out.ogr',flags = "overwrite", format='PostgreSQL', input=paste("Streams_WR", wrid, "_network", sep=''), type='line', output=paste("PG:dbname=",pg_db,sep=""), output_layer=paste("temporary.\"Streams_WR", wrid, "_network\"", sep=''), lco=c("FID=cat","SRID=EPSG:25833","GEOMETRY_NAME=geom"))



