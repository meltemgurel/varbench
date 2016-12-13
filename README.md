# varbench
varbench is a ....
### running varbench
clone workflow into working directory
```sh
$ git clone https://github.com/meltemgurel/varbench.git path/to/workdir
$ cd path/to/workdir
```
edit config specifying AT LEAST the sample file
```sh
$ vim config.yml
```
install dependencies into isolated environment
```sh
$ conda create -n myworkflow --file environment.yml
```
activate environment
```sh
$ source activate myworkflow
```
execute workflow
```sh
$ snakemake --configfile config.yml
```
