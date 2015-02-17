import os
import subprocess
import datetime
from optparse import OptionParser

parser = OptionParser()
parser.usage = '''
'''
parser.add_option("-Q"  , "--queue"  , dest = "queue"  , help = "choose queue. Default is shortcms"           , default = 'shortcms'                    )
parser.add_option("-N"  , "--njobs"  , dest = "njobs"  , help = "choose number of jobs.Default is -1"         , default = -1                            )
parser.add_option("-F"  , "--nfiles" , dest = "nfiles" , help = "choose number of files per job.Default is 1" , default =  1                            )
parser.add_option("-S"  , "--start"  , dest = "start"  , help = "choose starting file"                        , default =  0                            )
parser.add_option("-i"  , "--input"  , dest = "input"  , help = "input file list"                             , default =  "fileList.txt"               )
parser.add_option("-c"  , "--cfg"    , dest = "cfg"    , help = "input cfg file"                              , default =  "namelistSkimTemplate.dat"   )
parser.add_option("-d"  , "--outdir" , dest = "outdir" , help = "output dir on gwteray"                       , default =  "SkimJpsiPi_2012C"           )

(options,args) = parser.parse_args()  

njobs  = options.njobs
nfiles = options.nfiles  

njobs  = int(njobs)
nfiles = int(nfiles)

myworkingfolder = os.getcwd()

myTime = []
for i in list(datetime.datetime.now().timetuple())[:5] :
  myTime.append(str(i))
newFolder = '/gwteray/users/fiorendi/data2012/cfgs_{TIME}_{OUT}'.format(TIME='_'.join(myTime), OUT=options.outdir)
inputlist = myworkingfolder + '/' + options.input
subprocess.check_call(['mkdir', newFolder])

flist   = open(inputlist)
infiles = flist.readlines()

myFirstFile = int(options.start)
if myFirstFile > 1:
  infiles = infiles[int(myFirstFile)-1:]

if int(njobs) == -1 : 
  njobs = len(infiles)/int(nfiles) + 1*(len(infiles)%int(nfiles)>0) 
  
for j in range(njobs): 
  k = myFirstFile + j
  print " k is " + str(k)
  f   = open(myworkingfolder + '/' + options.cfg)
  f1  = open(newFolder + '/{M}_{NUM}.dat'.format(M=options.cfg.replace(".dat","").replace("../","").rstrip(), NUM=str(k)), "w")
  for line in f:
      newline    = None
      newnewline = None
      newfile    = None
      if 'inFileNameSub' in line:
          for i in infiles[k*nfiles:(k+1)*nfiles]:
            if i=="": continue
            newfile = line.replace('inFileNameSub', '').rstrip()
            print >> f1,  newfile.rstrip() + '  ' + i.rstrip()
            sample  = i.split("/")[7].rstrip() 
      if 'outputFileNameSub' in line:
          newline = line.replace('outputFileNameSub', '{TIME}/Analysis_{TAG}'.format(TIME=newFolder, NUM=str(k), TAG=sample )).rstrip()
      if newline:
        print >> f1,newline.strip()
      if 'outputTreeNameSub' in line:
          newnewline = line.replace('outputTreeNameSub', '{TIME}/Tree_{TAG}'.format(TIME=newFolder, NUM=str(k), TAG=sample )).rstrip()
      if newnewline:
        print >> f1,newnewline.strip()
      if not (newnewline or newline or newfile):
        print >> f1,line.rstrip() 

  f.close()
  f1.close()
  
  sh   = open("script.sh")
  shName = myworkingfolder + "/myScripts/script_{NUM}.sh".format(NUM=str(k))
  sh1  = open(shName,"w")
  subprocess.call(['chmod', '+x', shName])
  for shline in sh:
      mycfg = None
      myres = None
      if 'file.dat'  in shline:
        shline = shline.replace('file.dat', newFolder + '/{M}_{NUM}.dat'.format(M=options.cfg.replace(".dat","").rstrip(), NUM=str(k))).rstrip()
      if shline:
        print >> sh1,shline.strip()
      else: 
        print >> sh1, shline.rstrip()
  sh.close()
  sh1.close()
       
  subprocess.call(['echo','qsub', '-q', options.queue, shName ])
  subprocess.call(['qsub', '-q', options.queue, shName ])

