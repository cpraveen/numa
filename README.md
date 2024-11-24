# Numerical Analysis

See the book in html format here

https://cpraveen.github.io/numa

## Install 

```shell
conda install jupyter-book sphinx-exercise ghp-import
```

## Build html and publish

```shell
git checkout git@github.com:cpraveen/numa
jb build numa/
cd numa
ghp-import -n -p -f _build/html
```
