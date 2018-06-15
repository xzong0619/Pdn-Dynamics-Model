import os 
import numpy as np



class ProcstatOut(object):
    
    '''
    Handles data from procstat_output.txt
    '''
    
    fname = 'procstat_output.txt'

    def __init__(self):
    
        self.Spacing = None
        self.t = None
        self.events = None
    
    def ReadOut(self, fldr):
        '''
        Read procstat_output.txt
        '''
        MaxLen = np.int(2e4)
        with open(os.path.join(fldr, self.fname), 'r') as txt:
            RawTxt = txt.readlines()

        if len(RawTxt) - 1 > MaxLen * 3:  # Procstat uses 3 lines per outputs
            Spacing = np.int(np.floor((len(RawTxt) - 1)/(MaxLen*3)))
            RawTxt2 = []
            for i in range(0, MaxLen):
                RawTxt2.append(RawTxt[i*Spacing*3+1])
                RawTxt2.append(RawTxt[i*Spacing*3+2])
                RawTxt2.append(RawTxt[i*Spacing*3+3])
        else:
            Spacing = 1
            RawTxt2 = RawTxt[1:]

        t = []
        events = []
        for i in range(0, np.int(len(RawTxt2)/3)):
            t.append(np.float(RawTxt2[i*3].split()[3]))
            eventsTemp = RawTxt2[i*3+2].split()[1:]
            for j in range(0, np.int(len(eventsTemp))):
                eventsTemp[j] = np.int(eventsTemp[j])
            events.append(eventsTemp)

        self.Spacing = Spacing
        self.t = np.asarray(t)
        self.events = np.asarray(events)

Parent_dir = os.getcwd()
fldr = os.path.join(Parent_dir, 'zacros_inputs')
x = ProcstatOut()
x.ReadOut(fldr)
