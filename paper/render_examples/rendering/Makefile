MAPPER := /home/patrick/repos/transitmap/build/transitmap

all:
	cat render_example.json | $(MAPPER) -o render_example1.svg --expand-fronts=0 --render-node-circles=0 --render-stations=0 --line-width=15 --line-spacing=10 --render-node-connections=0
	cat render_example.json | $(MAPPER) -o render_example2.svg --expand-fronts=1 --render-node-fronts=1 --render-node-circles=0 --render-stations=0 --line-width=15 --line-spacing=10 --render-node-connections=0
	cat render_example.json | $(MAPPER) -o render_example3.svg --expand-fronts=1 --render-node-fronts=1 --render-node-circles=0 --render-stations=0 --line-width=15 --line-spacing=10 --render-node-connections=1
	cat render_example.json | $(MAPPER) -o render_example4.svg --expand-fronts=1 --render-node-circles=0 --render-stations=1 --line-width=15 --line-spacing=10 --render-node-connections=1
	inkscape render_example1.svg --export-pdf=render_example1.pdf
	inkscape render_example2.svg --export-pdf=render_example2.pdf
	inkscape render_example3.svg --export-pdf=render_example3.pdf
	inkscape render_example4.svg --export-pdf=render_example4.pdf
