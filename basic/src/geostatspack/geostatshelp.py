# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'geostatshelp.ui'
#
# Created: Thu Nov 13 11:18:05 2008
#      by: The PyQt User Interface Compiler (pyuic) 3.17.4
#
# WARNING! All changes made in this file will be lost!


from qt import *
##############################################################################
# geostatshelp.py
# Part of PAth Search Tool bAsed on Flexible Atomistic Reaction Image ANalysis
# (c) 2005-2008 by Jan M. Knaup , Knaup@bccms.uni-bremen.de
# all rights reserved
##############################################################################
# Licensed under the Non-Profit Open Software License version 3.0
# see file LICENSE for details.
##############################################################################

image0_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x14\x00\x00\x00\x14" \
    "\x08\x06\x00\x00\x00\x8d\x89\x1d\x0d\x00\x00\x02" \
    "\x0e\x49\x44\x41\x54\x38\x8d\xad\x94\xbf\x4b\x5b" \
    "\x51\x14\xc7\x3f\x37\xbe\x67\x22\xf1\xbd\xa4\xa5" \
    "\x04\x89\x10\x77\x71\x2e\xe4\x2f\x70\x34\xb8\x5a" \
    "\x87\x2e\x2e\xfe\x01\x0e\x6e\xce\x4a\xa6\xd0\x35" \
    "\xae\x01\xc1\xb5\xb4\x5d\xb2\x26\x10\x44\x0d\x91" \
    "\x2c\x31\xf5\xa5\xf2\xe8\x7b\x49\xa3\xc1\x34\x2f" \
    "\xf6\x74\x89\xf9\x61\x94\xda\xd7\x7e\xe1\xc2\xbd" \
    "\xe7\xc2\xe7\x9c\xfb\x3d\x87\xab\x44\x04\x80\x7a" \
    "\xbd\x2e\x37\x37\x37\x18\x86\x41\x2c\x16\x23\x14" \
    "\x0a\x29\x7c\x28\x00\xa0\x94\x92\xe3\xe3\x63\xba" \
    "\xdd\x2e\xba\xae\x53\x2c\x16\x91\x87\x4c\x7e\x80" \
    "\x00\xb6\x6d\x63\xdb\x36\xa5\x52\x89\xd3\xd3\x53" \
    "\x9a\xcd\xa6\x1f\x1e\xda\xc3\xc6\xb2\x2c\xf2\xf9" \
    "\x3c\xbd\x5e\x0f\xd7\x75\xd9\xda\xda\xf2\x05\x0c" \
    "\x28\xa5\x04\x20\x95\x4a\xe1\x38\x0e\xb6\x6d\xb3" \
    "\xb1\xb1\xc1\xcc\xcc\x8c\x2f\x0f\x95\x88\xa0\x94" \
    "\x1a\x5a\x56\x2e\x97\xa9\x56\xab\x24\x93\x49\x62" \
    "\xb1\xd8\x5f\x43\x03\x8f\x03\xcb\xcb\xcb\x74\x2e" \
    "\x2f\xb1\x36\x37\xe9\xd4\xeb\xd2\x3e\x3f\x97\xfb" \
    "\x6e\xf7\xe5\x0d\x1a\x54\x26\xe3\xba\xaa\x54\xe4" \
    "\x83\xae\xcb\x67\x90\x33\xa5\xe4\x2c\x1a\x95\x7a" \
    "\x26\x23\x22\xc2\x9f\x96\xf6\x54\x92\x8f\xfb\xfb" \
    "\x04\x3d\x8f\x30\x30\xaf\x69\xbc\x0a\x87\x29\x6f" \
    "\x6f\xf3\x6b\x71\x51\x96\xd6\xd6\x86\x36\x34\x1a" \
    "\x0d\xb9\xb8\xb8\x60\x61\x61\x81\x50\x28\x44\x22" \
    "\x91\x18\x01\x4d\xd7\x04\xe0\x3a\xf8\x8d\x2f\x87" \
    "\x87\xbc\x05\xda\x9a\xc6\xfc\xec\x2c\xdf\x2d\x8b" \
    "\x0e\x50\xd8\xdb\xc3\xbc\xbb\x13\x56\x57\x01\xc8" \
    "\x66\xb3\xac\xac\xac\x10\x8f\xc7\x69\xb5\x5a\xdc" \
    "\xde\xde\x4e\x57\xd8\x6e\x36\xb9\xf7\x3c\x1c\xe0" \
    "\xab\xa6\xc1\xdc\x1c\xcd\x4e\x87\x36\x50\x3d\x39" \
    "\xc1\xda\xd9\x41\xae\xaf\x01\xa8\xd5\x6a\x04\x83" \
    "\x41\xfa\xfd\x3e\xae\xeb\x12\x89\x44\x50\x80\x00" \
    "\x18\x8e\x01\xc0\x8f\x68\x8b\x77\x4b\x4b\xbc\xb9" \
    "\xba\x22\x32\xe8\xda\x4f\x20\x33\xb8\x6f\xbf\x6e" \
    "\x0f\x9f\x9c\x4e\xa7\xa5\x52\xa9\xa0\xeb\x3a\x96" \
    "\x65\xb1\xbb\xbb\x3b\x6a\x8a\xe1\x18\x62\x38\x86" \
    "\x88\x88\x7c\x3a\x3a\x1a\x9e\xdf\x83\xac\x4f\xde" \
    "\x0f\x1b\xe0\x79\x9e\xe4\x72\x39\x39\x38\x38\x90" \
    "\x42\xa1\x20\x22\x32\x9a\x43\x63\x54\xc1\x84\xa7" \
    "\x53\x96\x8c\x55\xf8\x94\xa6\x80\x2f\xd5\x73\xe0" \
    "\xa9\xc1\x7e\xa9\x4c\xd7\x14\xd3\x35\xa7\x06\xde" \
    "\x37\xf0\x39\xf0\x3f\x03\x1f\x6b\x38\x36\xff\x43" \
    "\x22\xa2\xc6\x7f\x9b\x09\x93\xfd\xc6\x7e\x03\x76" \
    "\x35\x47\x87\x2d\xbf\x85\xbe\x00\x00\x00\x00\x49" \
    "\x45\x4e\x44\xae\x42\x60\x82"

class GeostatsHelp(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)

        self.image0 = QPixmap()
        self.image0.loadFromData(image0_data,"PNG")
        if not name:
            self.setName("GeostatsHelp")

        self.setIcon(self.image0)

        GeostatsHelpLayout = QHBoxLayout(self,11,6,"GeostatsHelpLayout")

        self.helpTextBrowser = QTextBrowser(self,"helpTextBrowser")
        GeostatsHelpLayout.addWidget(self.helpTextBrowser)

        self.languageChange()

        self.resize(QSize(486,480).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)


    def languageChange(self):
        self.setCaption(self.__tr("Geostats Help"))
        self.setIconText(self.__tr("Help on Geostats"))
        QToolTip.add(self,self.__tr("Help on geostats.py"))
        QWhatsThis.add(self,self.__tr("Help on geostats.py"))
        self.helpTextBrowser.setText(self.__tr("<html>\n"
"<head>\n"
"<title>Geostats.py Help</title>\n"
"</head>\n"
"<body>\n"
"<h1>Geostats.py Help</h1>\n"
"<h2>Contents</h2>\n"
"<a href=\"#sec_intro\">Introduction</a><br>\n"
"<a href=\"#sec_formats\">File Formats</a><br>\n"
"<a href=\"#sec_menus\">Menus</a><br>\n"
"<a href=\"#sec_dialogs\">Dialogs</a><br>\n"
"<a href=\"#sec_copyright\">Copyright</a><br>\n"
"<a name=\"sec_intro\"><h2>Introduction</h2></a><br>\n"
"<p>Geostats.py is intended to to allow statistical analyses and simple geometry manipulations on atomistic models of molecules and solids.</p>\n"
"<p>Bug reports, feature suggestions and comments are welcome and should be addressed to the <a href=\"#sec_copyright\">author</a>.</p>\n"
"<a name=\"sec_formats\"><h2>Supported File Formats</h2></a>\n"
"<h3>input</h3>\n"
"<a name=\"fmt_gen\"><p><tt>.gen</tt>: Generic file format from old dftb implementations. Refer to the dftb documentation for reference.</p></a>\n"
"<a name=\"fmt_fmg\"><p><tt>.fmg</tt>: Flexible molecular geometry. XML format for molecular geometries and trajectories. More flexible and easily extensible but not so human readable.</p></a>\n"
"<a name=\"fmt_xyz\"><p><tt>.xyz</tt>: Plain and general purpose element symbol, X, Y, Z file format. Does not support supercells or additional atomic properties.</p></a>\n"
"<h3>output</h3>\n"
"<p><tt>.gen</tt>: <a href=\"#fmt_gen\">see above</a></p>\n"
"<p><tt>.fmg</tt>: <a href=\"#fmt_fmg\">see above</a></p>\n"
"<p><tt>.fdf</tt>: SIESTA flexible document format. Suitable to be included as geometry information into a SIESTA input file.</p>\n"
"<p><tt>.pdb</tt>: Protein DataBase format. Very limited support, all atoms are stored as HETATM records, segments correspond to geometry layers.</p>\n"
"<p><tt>.xyz</tt>: <a href=\"#fmt_xyz\">see above.</a></p>\n"
"<p><tt>.xyzq</tt>: Very simple X, Y, Z Charge file format to specify point charges.</p>\n"
"<a name=\"sec_menus\"><h2>Menus</h2></a>\n"
"<h3>File</h3>\n"
"File operations\n"
"<h4></h4>\n"
"<h4>open</h4>\n"
"<p>Calls an open file dialog to read a geometry file for analysis or editing</p>.\n"
"<h4>save</h4>\n"
"<p>Saves to the last opened filename.</p>\n"
"<h4>save as</h4>\n"
"<p>Calls a save file dialog to choose filename, directory and filetype to save the current geometry.</p>\n"
"<h4>view in VMD</h4>\n"
"<p><i>Not yet implemented</i>Opens an external VMD instance to view or edit the current Geometry.</p>\n"
"<h3>Edit</h3>\n"
"Geometry editing operations\n"
"<h4>edit atoms</h4>\n"
"<p>Opens the edit atoms dialog, which allows to modify layers and edit all atomic properties known to the .fmg file format.</p>\n"
"<h4>periodic expand</h4>\n"
"<p><i>only available for supercell geometries.</i> Opens the periodic expand dialog which allows to expand the geometry along the supercell vectors.</p>\n"
"<h4>generate linkatoms</h4>\n"
"<p>Adds link-hydrogens to the DEFAULT layer of the geometry along all bonds between the DEFAULT and PCHR layers.</p>\n"
"<h3>Statistics</h3>\n"
"<h4>coordinations</h4>\n"
"<p>Adds atomic cooridination statistics to the main text area.</p>\n"
"<h4>charges</h4>\n"
"<p>Adds atomic charge statistics to the main text area.</p>\n"
"<h4>graphs</h4>\n"
"<p>Submenu from which wizards for graphing statistical functions can be opened.</p>\n"
"<h5>available wizards</h5>\n"
"<ul>\n"
"<li><b>element charge histograms</b> Creates one separate atomic charge histogram for each element.</li>\n"
"<li><b>radial distribution function</b> Calculates the radial distribution function for the geometry.</li>\n"
"</ul>\n"
"<h4>bond list</h4>\n"
"<p>Adds a complete bond table to the main text area.<i>Can generate large amounts of text.</i></p>\n"
"<h4>save statistics</h4>\n"
"<p>Saves the contents of the main text area to geostats.html in the directory from which geostats.py was called.</p>\n"
"<h3>Help</h3>\n"
"<h4>contents</h4>\n"
"<p>This dialog</p>\n"
"<h4>about</h4>\n"
"<p>A very shirt description of geostats.py</p>\n"
"<a name=\"sec_dialogs\"><h2>Dialogs</h2>\n"
"<h3>edit atoms</h3>\n"
"<p>The edit atoms provides two functions: 1. a table view of the atoms, in which position element, position, layer, charge and subtype of each atom can be edited.</p>\n"
"<p>Additionally, atoms can deleted from and added to the tables, using buttons at the lower right corner of the table.</p>\n"
"<p>2. a table of geometry layers, in which existing layers can be deleted and new layers can be created. If a layer is deleted, all atoms in that layer are automatically assigned to the default layer.</p>\n"
"<p></p>\n"
"<h3>periodic expand</h3>\n"
"<p>The dialog displays the a,b,c supercell vectors in angstroms and allows to choose the number of old supercells to be in included in the new supercell in each direction. The old supercell itself would be 1,1,1. To double the supercell in all directions, one would enter 2,2,2.</p>\n"
"<p>Periodic expansion is inly available,  if the loaded geometry is specified as periodic and has supercell vectors.</p>\n"
"<h3>element charges histograms</h3>\n"
"<p>The element charges histograms wizard applies a binning procedure to the atomic charges of all atoms present int he geometry. The histograms are separate for each element.</p>\n"
"<p>In step 1, choose the number of bins per histogram. Do not choose a too fine binning grid.</p>\n"
"<p>In step 2, enter the filename into which the gnuplot-friendly histogram data should be written. Also, if gnuplot and the Gnuplot python module are available on your system, choose wheter to call gnupot directly to display the histograms and optionally save the graph in one of the following file formats:\n"
"<ul>\n"
"<li><tt>.eps</tt>: encapsulated postscript</li>\n"
"<li><tt>.png</tt>: portable network graphics</li>\n"
"<li><tt>.emf</tt>: enhanced metafile</li>\n"
"<li><tt>.fig</tt>: xfig graphics file</li>\n"
"</ul>\n"
"</p>\n"
"<p>in step 3, just klick finish to start the generation of the historgams. Progress bars inform about the progress of the binning progress over the lements.</p>\n"
"<h3>radial distribution function</h3>\n"
"<p>The radial distribution function wizard calculates the rdf of all atoms in the loaded geometry. for Molecules, the cutoff radius and density are derived from a sphere, tightly circumscribing the molecule. For a supercell structure these data are derived from the cell vectors.</p>\n"
"<p>In step 1, the enter stepwidth of the rdf calculation.</p>\n"
"<p>In step 2, enter the filename into which the gnuplot-friendly rdf  data should be written. Also, if gnuplot and the Gnuplot python module are available on your system, choose wheter to call gnupot directly to display the rdf and optionally save the graph in one of the following file formats:\n"
"<ul>\n"
"<li><tt>.eps</tt>: encapsulated postscript</li>\n"
"<li><tt>.png</tt>: portable network graphics</li>\n"
"<li><tt>.emf</tt>: enhanced metafile</li>\n"
"<li><tt>.fig</tt>: xfig graphics file</li>\n"
"</ul>\n"
"</p>\n"
"<p>In step 3, just klick finish to start the calculation of the rdf. Progress bars inform about the progress of the bond length binning and normalization.</p>\n"
"<p>Normalization is shaky and will not work as expected for molecules or surface slab models.</p>\n"
"<a name=\"sec_copyright\"><h2>Copyright</h2></a>\n"
"<p>This program was written by <a href=\"mailto:Jan.Knaup@bccms.uni-bremen.de\"><b>Jan M. Knaup</b> &lt;Jan.Knaup@bccms.uni-bremen.de&gt;</a> of the <a href=\"http://www.bccms.uni-bremen.de\"><b>Bremen Center for Computational Materials Science</b></a>. It may be redistributed without charge, as long as this copyright notice is left intact. Any software derived from this package must contain a reference to geostats.py and proper mention of the copyrights on geostats.py</p>\n"
"</body>"))
        QToolTip.add(self.helpTextBrowser,self.__tr("Help on geostats.py"))
        QWhatsThis.add(self.helpTextBrowser,self.__tr("Help on geostats.py"))


    def __tr(self,s,c = None):
        return qApp.translate("GeostatsHelp",s,c)
