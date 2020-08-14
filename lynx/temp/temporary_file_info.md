Temp files ending with `.csv`/`.xlsx` will be checked every time when LipidLynxX api server starts.

**if**

some temp files that older than defined lifetime threshold

**or**

if there are more files than the `temp_max_files` value

Files fit to **any of above situation** will be removed.

Examples:

+ If 200 files were generated in 2 days while 
`temp_max_days` set to `3` days and 
`temp_max_files` is set to `99` files,
only latest 99 files will be kept and all other older files will be removed.

+ If 40 files were generated in 2 days while 
`temp_max_days` set to `3` days and 
`temp_max_files` is set to `99` files,
all 40 files will be kept.


Please modify the settings in `LipidLynxX/lynx/config.ini`

Change following settings:

```
temp_folder = lynx/temp
temp_max_days = 3
temp_max_files = 99
```


`temp_folder` is set to this folder by default. 
Any other folder **MUST be provided as absolute path**.

`temp_max_days` is set to `3` days by default. 
Accept `int` values only.

`temp_max_files` is set to `99` by default. 
Accept `int` values only.
