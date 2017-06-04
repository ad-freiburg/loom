MAPPER := /home/patrick/repos/transitmap/build/transitmap

all:
	cat periphery_example1.json | $(MAPPER) -o periphery_example1.svg --splitting-optim=1 --simple-station-render-heur=0 --collapse-line-partners=1 --render-node-fronts=0 --expand-fronts=1 --render-node-circles=0 --render-stations=1 --line-width=15 --line-spacing=10 --render-node-connections=1
	inkscape periphery_example1.svg --export-pdf=periphery_example1.pdf

	cat periphery_example2.json | $(MAPPER) -o periphery_example2.svg --splitting-optim=1 --simple-station-render-heur=0 --collapse-line-partners=1 --render-node-fronts=0 --expand-fronts=1 --render-node-circles=0 --render-stations=1 --line-width=15 --line-spacing=10 --render-node-connections=1
	inkscape periphery_example2.svg --export-pdf=periphery_example2.pdf
