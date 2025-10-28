SRC=$(wildcard *.md)
PDF=$(SRC:.md=.pdf)

all: $(PDF)

%.pdf: %.md
	myst build $< --execute --pdf

pdf:
	myst build --execute --pdf

clean:
	myst clean
