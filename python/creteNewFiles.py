import os
from parallelize import parallelize

inputpath='root://maite.iihe.ac.be:1094//pnfs/iihe/cms/store/user/anmalara/Inputs'
outputhpath=os.environ['LEAFPATH']+'/StudyGroup/Inputs'
extra='standard=True extrajets=True pfcands=False triggerobjects=False stablegenparticles=False'
type='DATA'
start='0'
stop='1000'
config=os.environ['ANALYZERPATH']+'/python/ntuplizer_cfg.py'
years = ['UL2018','UL2017','UL2016','EOY2018D','EOY2018ABC','EOY2017','EOY2016']
srs = ['ER1']
cats = ['tautau','mutau','etau']

commands = []
for sr in srs:
    for cat in cats:
        for year in years:
            year2 = year.replace('20','').replace('EOY','20').replace('D','').replace('ABC','')
            fname = cat+'_0b_'+sr+'_'+year
            cmd = 'cmsRun '+config+' type='+type+' infilename='+inputpath+'/MINIAOD_'+fname+'.root outfilename='+outputhpath+'/NTuples_'+fname+'.root idxStart='+start+' idxStop='+stop+' year='+year2+' '+extra
            # print(cmd)
            commands.append([cmd])

a = parallelize(commands, ncores=4, remove_temp_files=True)
