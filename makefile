all: pdf

pdf:
	myst build --execute --pdf

clean:
	myst clean
