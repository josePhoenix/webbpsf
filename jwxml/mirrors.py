from .utils import iterchildren, etree

import logging
_log = logging.getLogger('jwxml')


#---------------------------------------------------------------------------------
#   Mirror Move related classes

class Segment_Update(object):
    """ Class for representing one single mirror update (will be inside of groups in SURs)
    """
    def __init__(self, xmlnode):
        if xmlnode.attrib['type'] != 'pose': raise NotImplemented("Only Pose updates supported yet")

        self.id = int(xmlnode.attrib['id'])
        self.type = xmlnode.attrib['type']
        self.segment = xmlnode.attrib['seg_id'][0:2]
        self.absolute = xmlnode.attrib['absolute'] =='true'
        self.coord= xmlnode.attrib['coord'] #local or global
        self.stage_type= xmlnode.attrib['stage_type']  # recenter_fine, fine_only, none

        self.units = dict()
        self.moves = dict()
        for move in iterchildren(xmlnode):
            #print(move.tag, move.text )
            self.moves[move.tag] =float(move.text)
            self.units[move.tag] = move.attrib['units']
            #X_TRANS, Y_TRANS, PISTON, X_TILT, Y_TILT, CLOCK
        #allowable units:
		#units="id"
		#units="meters"
		#units="none"
		#units="radians"
		#units="sag"
		#units="steps"
		#
        # pose moves will only ever have meters/radians as units
    def __str__(self):
        return ("Update %d, %s, %s: "% (self.id, 'absolute' if self.absolute else 'relative', self.coord)) + str(self.moves)
    def shortstr(self):
        outstr = ("Update %d: %s, %s, %s {"% (self.id, self.segment, 'absolute' if self.absolute else 'relative', self.coord))

        outstr+= ", ".join([ coordname+"=%.3g" % self.moves[coordname] for coordname in ['PISTON','X_TRANS','Y_TRANS','CLOCK', 'X_TILT','Y_TILT']])
        #for coordname in ['PISTON','X_TRANS','Y_TRANS','CLOCK', 'X_TILT','Y_TILT']:
            #outstr+=coordname+"=%.3g" % self.moves[coordname]
        outstr+="}"
        return outstr

    @property
    def xmltext(self):
        """ The XML text representation of a given move """
        text= '        <UPDATE id="{0.id}" type="{0.type}" seg_id="{0.segment}" absolute="{absolute}" coord="{0.coord}" stage_type="{0.stage_type}">\n'.format( self, absolute = str(self.absolute).lower())
        for key in ['X_TRANS','Y_TRANS','PISTON','X_TILT', 'Y_TILT', 'CLOCK']:
            if key in self.moves:
                text+='            <{key}  units="{unit}">{val:E}</{key}>\n'.format(key=key, unit=self.units[key], val=self.moves[key])
        text+= '        </UPDATE>\n'
        return text

    def toGlobal(self):
        """ Return moves cast to global coordinates """
        if self.coord =='global':
            return self.moves
        else:
            raise NotImplemented("Error")


    def toLocal(self):
        """ Return moves cast to local coordinates """
        if self.coord =='local':
            return self.moves
        else:
            raise NotImplemented("Error")
            # TO implement based on Ball's 'pmglobal_to_seg' in ./wfsc_core_algs/was_core_pmglobal_to_seg.pro
            # or the code in ./segment_control/mcs_hexapod_obj__define.pro


class SUR(object):
    """ Class for parsing/manipulating Segment Update Request files

    """
    def __init__(self, filename):
        """ Read a SUR from disk """
        self.filename=filename

        self._tree = etree.parse(filename)

        for tag in ['creator','date','time','version', 'operational']:
            self.__dict__[tag] = self._tree.getroot().attrib[tag]
        for element in self._tree.getroot().iter():
            if element.tag =='CONFIGURATION_NAME':  self.configuration_name = element.text
            if element.tag =='CORRECTION_ID':  self.correction_id = element.text

        self.groups = []
        for grp in self._tree.getroot().iter('GROUP'):
            myupdates = []
            for update in grp.iter('UPDATE'):
                myupdates.append(Segment_Update(update))
            self.groups.append(myupdates)

    def __str__(self):
        outstr = "SUR %s\n" % self.filename #, type=%s, coords=%s\n" % (self.filename, 'absolute' if self.absolute else 'relative', self.coord)
        for igrp, grp in enumerate(self.groups):
            outstr+= "\tGroup %d\n" % (igrp+1)
            for update in grp:
                outstr+= "\t\t"+str(update)+"\n"
        return outstr

    @property
    def xmltext(self):
        """ The XML text representation of a given move """
        text = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<SEGMENT_UPDATE_REQUEST creator="?" date="{date}" time="{time}" version="0.0.1" operational="false" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="../../setup_files/schema/segment_update_request.xsd">
    <CONFIGURATION_NAME>{self.configuration_name}</CONFIGURATION_NAME>
    <CORRECTION_ID>{self.correction_id}</CORRECTION_ID>\n""".format(self=self, date='YYYY-MM-DD', time='HH:MM:SS')
    # FIXME add date and time keywords for real
        for igrp, grp in enumerate(self.groups):
            text+='    <GROUP id="{id}">\n'.format(id=igrp+1)
            for update in grp:
                text+=update.xmltext
            text+='    </GROUP>\n'
        text+= '</SEGMENT_UPDATE_REQUEST>'
        return text


    #@property
    #def name(self): return self._tree.getroot().attrib['name']
