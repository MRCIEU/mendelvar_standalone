## Initial set-up for standalone MendelVar usage.
```
mendelvar_intial_setup.md
```
Manual changes required (see the "Populate OMIM data" section in the `mendelvar_intial_setup.md` file)

Default location to clone the repository with MendelVar scripts into: `$HOME/bin/mendelvar_standalone`

**Make sure that two default locations below exist:**

Default location for MendelVar data: `$HOME/MendelVar`

Default location for MendelVar results: `$HOME/MendelVar_out`

## Script with MendelVar database updates. 
```
sh ./user_input/2_mendelvar_regular_update.sh
```
Manual changes required (see lines 13, 30, 33 in the script)
## Script for runs of standalone MendelVar jobs. 
```
sh ./user_input/1_run_user_input.sh
```
Manual changes required (see line 15 in the script)

## MendelVar tutorial
https://www.notion.so/MendelVar-tutorial-ab91d2a6acb846f2b9f2978fcd942dd5
