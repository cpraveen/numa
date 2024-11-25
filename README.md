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

## Docs

* [MyST syntax cheat sheet](https://jupyterbook.org/en/stable/reference/cheatsheet.html)
* [Proof, Theorems, Algorithms](https://jupyterbook.org/en/stable/content/proof.html)
* [Exercise](https://mystmd.org/guide/exercises)
