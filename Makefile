julia-shell:
	docker exec -it csc530-project_julia_1 /bin/bash

generate:
	julia generate/call.jl --t=10 --l=100 --k=8 --d=1 --gc=50 --m=2 --output-dir=../data/test

pairmotif:
	julia pairmotif/call.jl --dna-file=../data/test/dna.json --k=8 --d-threshold=1 --output-dir=../data/test

pairmotif-fast:
	julia pairmotif-faster/call.jl --dna-file=../data/test/dna.json --k=8 --n=60 --use-prob=false --d-threshold=1 --output-dir=../data/test

evaluate:
	julia evaluate.jl --eval_file=../data/eval_pairmotif.json --temp-dir=../data/test