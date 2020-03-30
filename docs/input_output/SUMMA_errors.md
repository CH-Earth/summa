# SUMMA Debugging

<a id="error_message_structure"></a>
## Error message structure 
SUMMA will prints error messages that can help diagnosing why a run has failed.  Error messages are concatenated strings showing SUMMA's progress through different functions until it reached the point of failure. An example error message might read:

```
summa_readForcing/read_force/getFirstTimestep/first requested simulation timestep not in any forcing file
```

The part before the indicates SUMMA's progress through different functions. The final item specifies why executing the final function failed. A straightforward way to find which routine created this message is to use `grep` (arguments: case `i`nsensitive, `r`ecursive search, line `n`umbers):

```
grep -irn "getFirstTimestep" ~/summa/
/home/wknoben/summa/build/source/engine/read_force.f90:122:  call getFirstTimestep(currentJulday,iFile,iRead,ncid,err,cmessage)
/home/wknoben/summa/build/source/engine/read_force.f90:233: subroutine getFirstTimestep(currentJulday,iFile,iRead,ncid,err,message)
/home/wknoben/summa/build/source/engine/read_force.f90:259: err=0; message="getFirstTimestep/"
/home/wknoben/summa/build/source/engine/read_force.f90:325: end subroutine getFirstTimestep
```
```
grep -irn "first requested simulation timestep not in any forcing file" ~/summa/
~/summa/build/source/engine/read_force.f90:319:   if(iFile==nFiles)then; err=99; message=trim(message)//'first requested simulation timestep not in any forcing file'; return; end if
```

The main routine is `read_force`, sub-routine is `getFirstTimestep` and the specific error message was generated on line 319 in `read_force.f90`.

<a id="common_errors"></a>
## Common errors

**_first requested simulation timestep not in any forcing file:_** First forcing time step does not match first simulation time step
```
summa_readForcing/read_force/getFirstTimestep/first requested simulation timestep not in any forcing file
```
*Solution:* Make sure that `simulStart` in `modelDecisions.txt` matches the first time step of the forcing data.


**Problem description**
```
summa_readForcing/readForcingData/problem reading forcing data: airpres/NetCDF: Start+count exceeds dimension boundgetFirstTimestep/
```
Solution