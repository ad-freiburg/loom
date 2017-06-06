TRANSITMAP_BIN := ../build/transitmap
OPTIM := ilp_impr
CORE_GRAPH := 1
SPLITTING_OPTIM := 1
CBC := 1
CBC_BIN := /home/patrick/repos/Cbc-2.9/build/bin/cbc

MAPS := maps/freiburg maps/chicago maps/stuttgart maps/turin maps/nyc_subway maps/busses_manhattan maps/dallas maps/boston

PARAMS := --coin-cbc-binary=$(CBC_BIN) --render-station-names=1 --glpk-time-limit=100000000 --output-stats=1  --render-node-fronts=1
RPARAMS := --coin-cbc-binary=$(CBC_BIN) --outline-width=1 --render-edges=0 --render-station-names=0 --glpk-time-limit=100000000 --output-stats=0  --render-node-fronts=0 --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH}


freiburg/network.svg: freiburg/network.json
	cat $< | $(TRANSITMAP_BIN) -o $@ $(PARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH} | sed '/stats/w freiburg/stats'

chicago/network.svg: chicago/network.json
	cat $< | $(TRANSITMAP_BIN) -o $@ $(PARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH} | sed '/stats/w chicago/stats'

dallas/network.svg: dallas/network.json
	cat $< | $(TRANSITMAP_BIN) -o $@ $(PARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH} | sed '/stats/w dallas/stats'

stuttgart/network.svg: stuttgart/network.json
	cat $< | $(TRANSITMAP_BIN) -o $@ $(PARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH} | sed '/stats/w stuttgart/stats'

nyc_subway/network.svg: nyc_subway/network.json
	cat $< | $(TRANSITMAP_BIN) -o $@ $(PARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH} | sed '/stats/w nyc_subway/stats'

busses_manhattan/network.svg: busses_manhattan/network.json
	cat $< | $(TRANSITMAP_BIN) -o $@ $(PARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH} | sed '/stats/w busses_manhattan/stats'

clean:
	rm -f freiburg/network.svg
	rm -f freiburg/stats
	rm -f chicago/network.svg
	rm -f chicago/stats
	rm -f dallas/network.svg
	rm -f dallas/stats
	rm -f stuttgart/network.svg
	rm -f stuttgart/stats
	rm -f nyc_subway/network.svg
	rm -f nyc_subway/stats
	rm -f busses_manhattan/network.svg
	rm -f busses_manhattan/stats

	rm -rf $(MAPS)

allmaps: $(MAPS)

maps/%: %/network.json
	mkdir -p $@
	mkdir -p $@/export/
	mkdir -p $@/z16/
	mkdir -p $@/z15/
	mkdir -p $@/z14/
	mkdir -p $@/z13/
	mkdir -p $@/z12/
	mkdir -p $@/z11/

	cp $< $@/export/

	# the reference used for download as3
	cat $< | $(TRANSITMAP_BIN) --coin-cbc-binary=$(CBC_BIN) --outline-width=3 --render-stats=1 --line-width=20 --line-spacing=10 -o $@/export/network.svg --glpk-mps-output-path=$@/export/network.mps --world-file-path=$@/export/  --resolution=0.1 --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --name=$(shell basename $@ | tr /a-z/ /A-Z/) --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH}
	inkscape $@/export/network.svg  --export-pdf=$@/export/network.pdf

	# z16
	# cat $< | $(TRANSITMAP_BIN) $(RPARAMS) --render-engine=svg_sep --line-width=10 --line-spacing=5 -o $@/z16/ --world-file-path=$@/z16/network_e.pngw  --resolution=0.5
	# cp $@/z16/network_e.pngw $@/z16/network_s.pngw
	# cp $@/z16/network_e.pngw $@/z16/network_n.pngw
	# inkscape $@/z16/edges.svg --export-png=$@/z16/network_e.png
	# inkscape $@/z16/nodes.svg --export-png=$@/z16/network_n.png
	# inkscape $@/z16/stations.svg --export-png=$@/z16/network_s.png
	# python scripts/gdal2tilesp.py -w none -z 16 -s EPSG:4326 $@/z16/network_e.png $@/edges
	# python scripts/gdal2tilesp.py -w none -z 16 -s EPSG:4326 $@/z16/network_s.png $@/stations
	# python scripts/gdal2tilesp.py -w none -z 16 -s EPSG:4326 $@/z16/network_n.png $@/nodes

	# z15
	cat $< | $(TRANSITMAP_BIN) $(RPARAMS) --render-engine=svg_sep --line-width=20 --line-spacing=10 -o $@/z15/ --world-file-path=$@/z15/network_e.pngw  --resolution=0.25
	cp $@/z15/network_e.pngw $@/z15/network_s.pngw
	cp $@/z15/network_e.pngw $@/z15/network_n.pngw
	inkscape $@/z15/edges.svg --export-png=$@/z15/network_e.png
	inkscape $@/z15/nodes.svg --export-png=$@/z15/network_n.png
	inkscape $@/z15/stations.svg --export-png=$@/z15/network_s.png
	python scripts/gdal2tilesp.py -w none -z 15 -s EPSG:4326 $@/z15/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 15 -s EPSG:4326 $@/z15/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 15 -s EPSG:4326 $@/z15/network_n.png $@/nodes


	# z14
	cat $< | $(TRANSITMAP_BIN) $(RPARAMS) --render-engine=svg_sep --line-width=40 --line-spacing=20 -o $@/z14/ --world-file-path=$@/z14/network_e.pngw  --resolution=0.125
	cp $@/z14/network_e.pngw $@/z14/network_s.pngw
	cp $@/z14/network_e.pngw $@/z14/network_n.pngw
	inkscape $@/z14/edges.svg --export-png=$@/z14/network_e.png
	inkscape $@/z14/nodes.svg --export-png=$@/z14/network_n.png
	inkscape $@/z14/stations.svg --export-png=$@/z14/network_s.png
	python scripts/gdal2tilesp.py -w none -z 14 -s EPSG:4326 $@/z14/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 14 -s EPSG:4326 $@/z14/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 14 -s EPSG:4326 $@/z14/network_n.png $@/nodes


	# z13
	cat $< | $(TRANSITMAP_BIN) $(RPARAMS) --render-engine=svg_sep --line-width=70 --line-spacing=35 -o $@/z13/ --world-file-path=$@/z13/network_e.pngw  --resolution=0.075
	cp $@/z13/network_e.pngw $@/z13/network_s.pngw
	cp $@/z13/network_e.pngw $@/z13/network_n.pngw
	inkscape $@/z13/edges.svg --export-png=$@/z13/network_e.png
	inkscape $@/z13/nodes.svg --export-png=$@/z13/network_n.png
	inkscape $@/z13/stations.svg --export-png=$@/z13/network_s.png
	python scripts/gdal2tilesp.py -w none -z 13 -s EPSG:4326 $@/z13/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 13 -s EPSG:4326 $@/z13/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 13 -s EPSG:4326 $@/z13/network_n.png $@/nodes


	# z12
	cat $< | $(TRANSITMAP_BIN) $(RPARAMS) --render-engine=svg_sep --line-width=90 --line-spacing=40 -o $@/z12/ --world-file-path=$@/z12/network_e.pngw  --resolution=0.03
	cp $@/z12/network_e.pngw $@/z12/network_s.pngw
	cp $@/z12/network_e.pngw $@/z12/network_n.pngw
	inkscape $@/z12/edges.svg --export-png=$@/z12/network_e.png
	inkscape $@/z12/nodes.svg --export-png=$@/z12/network_n.png
	inkscape $@/z12/stations.svg --export-png=$@/z12/network_s.png
	python scripts/gdal2tilesp.py -w none -z 12 -s EPSG:4326 $@/z12/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 12 -s EPSG:4326 $@/z12/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 12 -s EPSG:4326 $@/z12/network_n.png $@/nodes


	# z11
	cat $< | $(TRANSITMAP_BIN) $(RPARAMS) --render-engine=svg_sep --line-width=170 --line-spacing=70 -o $@/z11/ --world-file-path=$@/z11/network_e.pngw  --resolution=0.015
	cp $@/z11/network_e.pngw $@/z11/network_s.pngw
	cp $@/z11/network_e.pngw $@/z11/network_n.pngw
	inkscape $@/z11/edges.svg --export-png=$@/z11/network_e.png
	inkscape $@/z11/nodes.svg --export-png=$@/z11/network_n.png
	inkscape $@/z11/stations.svg --export-png=$@/z11/network_s.png

	# something is not right with multi-zoomlevel tile generationin parallel version, create them separately
	python scripts/gdal2tilesp.py -w none -z 11 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 11 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 11 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 10 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 10 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 10 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 9 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 9 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 9 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 8 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 8 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 8 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 7 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 7 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 7 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 6 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 6 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 6 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 5 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 5 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 5 -s EPSG:4326 $@/z11/network_n.png $@/nodes

	python scripts/gdal2tilesp.py -w none -z 4 -s EPSG:4326 $@/z11/network_e.png $@/edges
	python scripts/gdal2tilesp.py -w none -z 4 -s EPSG:4326 $@/z11/network_s.png $@/stations
	python scripts/gdal2tilesp.py -w none -z 4 -s EPSG:4326 $@/z11/network_n.png $@/nodes