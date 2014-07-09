#     +------------------------------------+
#     | Author: D. Menasce                 |
#     | Command line arguments parser      |
#     +------------------------------------+

import FWCore.ParameterSet.VarParsing as VarParsing
import inspect

#---------------------------------------------------------------------------------
def lineno():
    """Returns the current line number and module in the program."""
    msg = str(inspect.currentframe().f_back.f_lineno) 	      + \
          "]\t["                                      	      + \
	  inspect.currentframe().f_back.f_globals["__name__"] + \
	  "]"
    diff = 70 - len(msg)
    for n in range(diff) :
	msg += " "
    return msg

#---------------------------------------------------------------------------------
maxEvents = -1      # Setup reasonable defaults
verbose   = 0
dataType  = 'data'
output    = 'dataGrid.root'
useJSON   = 0
#---------------------------------------------------------------------------------
# Try parsing: use try-catch to trap abnd avoid execution during compilation stage 
# when a python file name has not yet been provided as input to cmsRun

try:
	options = VarParsing.VarParsing('standard')

        # Register two customized entries with the standard table
	options.register(
	                 'dataType',
			 'data',
			 VarParsing.VarParsing.multiplicity.singleton,
			 VarParsing.VarParsing.varType.string,
			 "Input data type (data or mc)"
	                )

	options.register(
	                 'useJSON',
			 0,
			 VarParsing.VarParsing.multiplicity.singleton,
			 VarParsing.VarParsing.varType.int,
			 "useJSON"
	                )
	
	options.register(
	                 'verbose',
			 0,
			 VarParsing.VarParsing.multiplicity.singleton,
			 VarParsing.VarParsing.varType.int,
			 "Verbosity"
	                )
	
	options.maxEvents = maxEvents 
	options.verbose   = verbose   
	options.dataType  = dataType  
	options.output    = output
	options.useJSON   = useJSON

        try:
		options.parseArguments()    # parse input arguments
	except:
		print lineno(), "]\tPlease ignore the above 'Failure' error: it is a fake"

	maxEvents = options.maxEvents
	verbose   = options.verbose  
	dataType  = options.dataType 
	output    = options.output
	useJSON   = options.useJSON
	useJSON   = 0
	
except:
	# Do nothing here: the purpose of catching this error is to avoid compilation errors
	print lineno(), " "

else:	

        print "\n"
	print "+------------------------------------------------------------+"
	print "|                                                            |"
	print "|                   Welcome to Bc2JpsiPi                     |"
	print "|                                                            |"
	print "+------------------------------------------------------------+"
        print "\n"
        
	
	print lineno(), "================ Processing " + \
	                str(maxEvents)                 + \
			" "                            + \
			dataType                       + \
			" events to "                  + \
			output                         + \
			" ========================"
	print "\n"
