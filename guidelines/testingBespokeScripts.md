## Unit test for bespoke scripts

These are temptative guidelines to write robust bespoke scripts. _Bespoke script_ are those
written for specific, one-off tasks. They are part of larger analyses and they are not finalized in
self contained packages. Bespoke scripts might be statististical analyses, plots, data formatting.
Formal unit tetsing for these scripts is usually not worth or easily feasible.

These guidelines use R as reference language.

* Use extensively `stopifnot(exp == obs)` statements. For example to test whether the number of rows
before and after an operation is as expected.

* Check Categories expressed as percentage data must sum up to 100

* Wrap code in functions and after a function definition apply a mini test. E.g.

```
cgPct<- function(x){
    ## Return %CG in string x
    cg_cnt<- nchar(gsub('A|T', '', toupper(x)))
    at_cnt<- nchar(gsub('C|G', '', toupper(x)))
    return(cg_cnt / (cg_cnt+at_cnt))
}
stopifnot(cgPct('atccggatat') == 0.4)
```

* Use extensively simulated data

* Select and identify columns in a table by **name** not by **position**. Rely on column position as little as possible.

### Managing data and files

* Keep backups but don't keep copies! A file should exist in one place only so if that file is changed or deleted that are no copies elsewhere out of sync.
Therefore prefer `mv` or `ln -s` over `cp`, likewise to copy between servers use `rsync --remove-source-files`.

* More files and analyses means more things to document and keep track of. Even if it doesn't cost anything, create files (and do analyses) only as needed.

* Mark temporary files and directories with `tmp` prefix or suffix. Delete these as soon as done.
