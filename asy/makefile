ASY = $(wildcard *.asy)
SVG = $(patsubst %.asy,%.svg,$(ASY))

all : $(SVG)

%.svg: %.asy
	asy -f svg $*.asy

clean:
	rm -f $(SVG)
