texinputs = TEXINPUTS=.:../../common:
bibinputs = BIBINPUTS=../../common:

procoptions = -pdf \
	 -output-directory=build \
	 -interaction=nonstopmode \
	 -halt-on-error

all: build/report.pdf

build/report.pdf: FORCE
	$(texinputs) $(bibinputs) latexmk $(procoptions) report.tex

preview:
	$(texinputs) $(bibinputs) latexmk $(procoptions) -pvc report.tex

FORCE:


build:
	mkdir -p build

clean:
	rm -rf build

.PHONY: all clean FORCE preview
