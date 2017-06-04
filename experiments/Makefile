TRANSITMAP_BIN := ../build/transitmap
OPTIM := ilp_impr
CORE_GRAPH := 1
SPLITTING_OPTIM := 1
CBC := 1

PARAMS := --render-station-names=1 --glpk-time-limit=100000000 --output-stats=1  --render-node-fronts=1
RPARAMS := --render-stations=1 --render-node-connections=1 --render-station-names=0 --glpk-time-limit=100000000 --output-stats=0  --render-node-fronts=0


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


maps/%: %/network.json
	mkdir $@
	mkdir $@/res10/
	# res 1/10
	cat $< | $(TRANSITMAP_BIN) -o $@/res10/tmp.svg --world-file-path=$@/res10.pngw  --resolution=0.1 $(RPARAMS) --use-cbc=${CBC} --splitting-optim=${SPLITTING_OPTIM} --optim-method=${OPTIM} --create-core-optim-graph=${CORE_GRAPH}
	inkscape $@/res10/tmp.svg --export-png=$@/res10.png --export-dpi=500