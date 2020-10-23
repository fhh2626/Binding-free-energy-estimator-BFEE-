##
## Binding Free Energy Estimator v1.1
##
## A plugin for automatical binding free energy calculation
##
## Email: fhh2626_at_gmail.com
## Date 2020.7.16
##
## cite: Haohao Fu, et al. JCIM, 2018, 58, 556
##
## Changelog:
##    2018.11.8 v0.5 beta:
##        H-mass repartitioning supported and numstep reduced
##        added water box based on the protein-ligand distance vector
##        reduced frequency of writing colvars outputs
##        now default pressure is 1.01325 bar
##    2018.11.12 v0.51 beta:
##        more strict check for polar_theta and polar_phi to guarantee the rotational invariance
##        by default, WTM-eABF instead classical eABF is used to accelerate sampling
##    2018.11.13 v0.52 beta:
##        minor big fixes
##    2018.11.20 v0.53 beta:
##        removed H-mass repartitioning, which is not stable in free-energy calculation
##        by default, colvarsRestartFrequency=outputFreq
##        reduced stepsPerCycle for stability
##    2020.6.12 v0.6 beta:
##        adapted the latest version of Colvars
##        uses stochastic velocity rescaling instead of Langevin barostat to improve efficiency
##    2020.6.15 v1.0:
##        alchemical free-energy calculation now supported
##    2020.6.16 v1.01:
##        minor bug fixes
##    2020.6.18 v1.02:
##        switched to bidirectional FEP by default
##    2020.7.16 v1.1:
##        now users can set reference and protein seperately
##        making tuning the pulling direction easily
##

package provide BFEEstimator 1.1

namespace eval BFEE {

variable version "v1.1"

variable w                                         ;# handle to main window

# below is for integration
# record the pmf
variable x []
variable G []
variable width 0
variable BOLTZMANN 0.0019872041

# for H-mass repartitioning
# hydrogen mass
variable HMass 1.0080
# distance threhold for covalent bonding
variable dThrehold 1.8

variable jobType 
variable temperature 
variable sel1
variable sel2
variable analyzermsd 
variable analyzeeulerTheta 
variable analyzeeulerPhi 
variable analyzeeulerPsi 
variable analyzepolarTheta
variable analyzepolarPhi 
variable analyzeTemperature 
variable r_star

variable psfDir
variable coorDir
variable velDir
variable xscDir
variable psfDirLig
variable coorDirLig
variable velDirLig
variable xscDirLig
variable paramDirList
variable rmsdDir
variable eulerThetaDir
variable eulerPhiDir
variable eulerPsiDir
variable polarThetaDir
variable polarPhiDir
variable rDir
variable unrmsdDir

# the GUI of BFEE plugin
proc ::BFEE::drawGUI {} {
	variable w
	variable jobType 
	variable temperature 
	variable sel1
	variable sel2
    variable sel3
	variable analyzermsd 
	variable analyzeeulerTheta 
	variable analyzeeulerPhi 
	variable analyzeeulerPsi 
	variable analyzepolarTheta
	variable analyzepolarPhi 
	variable analyzeTemperature 
	variable r_star
	variable version

	# If already initialized, just turn on
	if { [winfo exists .textview] } {
		wm deiconify $w
		return
	}

	set w [toplevel ".bfee"]
	wm title $w "Binding Free Energy Calculation"
	wm resizable $w 0 0


	ttk::label $w.discription -text "                  Binding Free Energy Calculation"  -borderwidth 5 -relief sunken -padding "5 5"

	ttk::radiobutton $w.jobType1 -text "Protein:Ligand" -variable ::BFEE::jobType -value proLig -width 20
	ttk::radiobutton $w.jobType2 -text "Protein:Protein" -variable ::BFEE::jobType -value proPro -state disabled
	set jobType proLig

	ttk::notebook $w.setupTab
	ttk::frame $w.setupTab.setup
	ttk::frame $w.setupTab.analyze

	# in setup notebook
	ttk::labelframe $w.setupTab.setup.input -text "Input for Complex"
	ttk::label $w.setupTab.setup.input.psfFile -text "psf file:"
	ttk::entry $w.setupTab.setup.input.psfEntry -textvariable ::BFEE::psfDir -width 21
	ttk::button $w.setupTab.setup.input.psfButton -text "Browse" -command ::BFEE::findPsfDir -width 9
	ttk::label $w.setupTab.setup.input.coorFile -text "coor file:"
	ttk::entry $w.setupTab.setup.input.coorEntry -textvariable ::BFEE::coorDir -width 21
	ttk::button $w.setupTab.setup.input.coorButton -text "Browse" -command ::BFEE::findCoorDir -width 9
	ttk::label $w.setupTab.setup.input.velFile -text "vel file:"
	ttk::entry $w.setupTab.setup.input.velEntry -textvariable ::BFEE::velDir -width 21
	ttk::button $w.setupTab.setup.input.velButton -text "Browse" -command ::BFEE::findVelDir -width 9
	ttk::label $w.setupTab.setup.input.xscFile -text "xsc file:"
	ttk::entry $w.setupTab.setup.input.xscEntry -textvariable ::BFEE::xscDir -width 21
	ttk::button $w.setupTab.setup.input.xscButton -text "Browse" -command ::BFEE::findXscDir -width 9

	#ttk::labelframe $w.setupTab.setup.inputLig -text "Input for Ligand"
	#ttk::label $w.setupTab.setup.inputLig.psfFile -text "psf file:"
	#ttk::entry $w.setupTab.setup.inputLig.psfEntry -textvariable ::BFEE::psfDirLig -width 21
	#ttk::button $w.setupTab.setup.inputLig.psfButton -text "Browse" -command ::BFEE::findPsfDirLig -width 9
	#ttk::label $w.setupTab.setup.inputLig.coorFile -text "coor file:"
	#ttk::entry $w.setupTab.setup.inputLig.coorEntry -textvariable ::BFEE::coorDirLig -width 21
	#ttk::button $w.setupTab.setup.inputLig.coorButton -text "Browse" -command ::BFEE::findCoorDirLig -width 9
	#ttk::label $w.setupTab.setup.inputLig.velFile -text "vel file:"
	#ttk::entry $w.setupTab.setup.inputLig.velEntry -textvariable ::BFEE::velDirLig -width 21
	#ttk::button $w.setupTab.setup.inputLig.velButton -text "Browse" -command ::BFEE::findVelDirLig -width 9
	#ttk::label $w.setupTab.setup.inputLig.xscFile -text "xsc file:"
	#ttk::entry $w.setupTab.setup.inputLig.xscEntry -textvariable ::BFEE::xscDirLig -width 21
	#ttk::button $w.setupTab.setup.inputLig.xscButton -text "Browse" -command ::BFEE::findXscDirLig -width 9

	ttk::labelframe $w.setupTab.setup.otherPar -text "Other parameters"
	ttk::label $w.setupTab.setup.otherPar.temp -text "Temperature:"
	ttk::entry $w.setupTab.setup.otherPar.tempEntry -textvariable ::BFEE::temperature
	set temperature 300
	ttk::label $w.setupTab.setup.otherPar.param -text "Parameter files:"
	listbox $w.setupTab.setup.otherPar.paramDir -listvariable ::BFEE::paramDirList -height 4 -width 21
	ttk::scrollbar $w.setupTab.setup.otherPar.scroll -command "$w.setupTab.setup.otherPar.paramDir xview" -orient horizontal
	$w.setupTab.setup.otherPar.paramDir configure -xscrollcommand "$w.setupTab.setup.otherPar.scroll set"
	ttk::button $w.setupTab.setup.otherPar.add -text "Add" -command ::BFEE::addParamDir -width 6
	ttk::button $w.setupTab.setup.otherPar.clear -text "Clear" -command ::BFEE::clearParamDir -width 6
    
    # reference is used to calculate Eular and spherical angles
    # and determine the distance and direction of pulling simulation
	ttk::label $w.setupTab.setup.otherPar.sel1 -text "Select reference:"
	ttk::entry $w.setupTab.setup.otherPar.sel1Entry -textvariable ::BFEE::sel1
	set sel1 "segname SH3D"
    # ligand should be the whole drug molecule
	ttk::label $w.setupTab.setup.otherPar.sel2 -text "Select ligand:"
	ttk::entry $w.setupTab.setup.otherPar.sel2Entry -textvariable ::BFEE::sel2
	set sel2 "segname PPRO"
    # protein should be the whole host molecule, used in step 8
    # (the protein will be removed)
    ttk::label $w.setupTab.setup.otherPar.sel3 -text "Select protein:"
	ttk::entry $w.setupTab.setup.otherPar.sel3Entry -textvariable ::BFEE::sel3
	set sel3 "segname SH3D"

	ttk::button $w.setupTab.setup.generate -text "Generate Inputs" -command ::BFEE::generateFiles

	ttk::labelframe $w.setupTab.setup.contactus -text "Readme"
	ttk::label $w.setupTab.setup.contactus.readme -text \
"step 1: Set up corresponding options and click
            <Generate Inputs>
step 2: Run simulations. Change <numstep> in 
            NAMD config file, <centers> of harmonic
            restraints and <lower-> and <upperboundary>
            of each collective variable if needed
step 3: Calculate binding free energies using 
            <Analyze> tab"
	#ttk::button $w.setupTab.setup.contactus.generate -text "About BFEEstimator" -command ::BFEE::aboutus
	#ttk::label $w.setupTab.setup.contactus.chipot -text "Chris Chipot: Christophe.Chipot@univ-lorraine.fr"
	#ttk::label $w.setupTab.setup.contactus.haohao -text "Haohao Fu: fhh2626@gmail.com"

	# in analyze notebook

	ttk::labelframe $w.setupTab.analyze.input -text "Input for PMFs (*.czar.pmf)"
	ttk::label $w.setupTab.analyze.input.boundstate -text "Bound state:"
	ttk::label $w.setupTab.analyze.input.rmsdFile -text "RMSD:" -anchor center
	ttk::entry $w.setupTab.analyze.input.rmsdEntry -textvariable ::BFEE::rmsdDir -width 23
	ttk::button $w.setupTab.analyze.input.rmsdButton -text "Browse" -command ::BFEE::findrmsdDir -width 9
	ttk::label $w.setupTab.analyze.input.eulerThetaFile -text "Theta:" -anchor center
	ttk::entry $w.setupTab.analyze.input.eulerThetaEntry -textvariable ::BFEE::eulerThetaDir -width 23
	ttk::button $w.setupTab.analyze.input.eulerThetaButton -text "Browse" -command ::BFEE::findeulerThetaDir -width 9
	ttk::label $w.setupTab.analyze.input.eulerPhiFile -text "Phi:" -anchor center
	ttk::entry $w.setupTab.analyze.input.eulerPhiEntry -textvariable ::BFEE::eulerPhiDir -width 23
	ttk::button $w.setupTab.analyze.input.eulerPhiButton -text "Browse" -command ::BFEE::findeulerPhiDir -width 9
	ttk::label $w.setupTab.analyze.input.eulerPsiFile -text "Psi:" -anchor center
	ttk::entry $w.setupTab.analyze.input.eulerPsiEntry -textvariable ::BFEE::eulerPsiDir -width 23
	ttk::button $w.setupTab.analyze.input.eulerPsiButton -text "Browse" -command ::BFEE::findeulerPsiDir -width 9
	ttk::label $w.setupTab.analyze.input.polarThetaFile -text "theta:" -anchor center
	ttk::entry $w.setupTab.analyze.input.polarThetaEntry -textvariable ::BFEE::polarThetaDir -width 23
	ttk::button $w.setupTab.analyze.input.polarThetaButton -text "Browse" -command ::BFEE::findpolarThetaDir -width 9
	ttk::label $w.setupTab.analyze.input.polarPhiFile -text "phi:" -anchor center
	ttk::entry $w.setupTab.analyze.input.polarPhiEntry -textvariable ::BFEE::polarPhiDir -width 23
	ttk::button $w.setupTab.analyze.input.polarPhiButton -text "Browse" -command ::BFEE::findpolarPhiDir -width 9
	ttk::label $w.setupTab.analyze.input.rFile -text "R:" -anchor center
	ttk::entry $w.setupTab.analyze.input.rEntry -textvariable ::BFEE::rDir -width 23
	ttk::button $w.setupTab.analyze.input.rButton -text "Browse" -command ::BFEE::findrDir -width 9

	ttk::label $w.setupTab.analyze.input.unboundstate -text "Unbound state:"
	ttk::label $w.setupTab.analyze.input.unrmsdFile -text "RMSD:" -anchor center
	ttk::entry $w.setupTab.analyze.input.unrmsdEntry -textvariable ::BFEE::unrmsdDir -width 23
	ttk::button $w.setupTab.analyze.input.unrmsdButton -text "Browse" -command ::BFEE::findunrmsdDir -width 9

	ttk::labelframe $w.setupTab.analyze.forceConst -text "Force Constants (in NAMD units)"
	ttk::label $w.setupTab.analyze.forceConst.bound -text "Bound state:"
	ttk::label $w.setupTab.analyze.forceConst.rmsd -text "RMSD:" -width 6 -anchor w
	ttk::entry $w.setupTab.analyze.forceConst.rmsdEntry -textvariable ::BFEE::analyzermsd -width 6
	set analyzermsd 10
	ttk::label $w.setupTab.analyze.forceConst.eulerTheta -text "Theta:" -width 6 -anchor center
	ttk::entry $w.setupTab.analyze.forceConst.eulerThetaEntry -textvariable ::BFEE::analyzeeulerTheta -width 6
	set analyzeeulerTheta 0.1
	ttk::label $w.setupTab.analyze.forceConst.eulerPhi -text "Phi:" -width 6 -anchor center
	ttk::entry $w.setupTab.analyze.forceConst.eulerPhiEntry -textvariable ::BFEE::analyzeeulerPhi -width 6
	set analyzeeulerPhi 0.1
	ttk::label $w.setupTab.analyze.forceConst.eulerPsi -text "Psi:" -width 6 -anchor w
	ttk::entry $w.setupTab.analyze.forceConst.eulerPsiEntry -textvariable ::BFEE::analyzeeulerPsi -width 6
	set analyzeeulerPsi 0.1
	ttk::label $w.setupTab.analyze.forceConst.polarTheta -text "theta:" -width 6 -anchor center
	ttk::entry $w.setupTab.analyze.forceConst.polarThetaEntry -textvariable ::BFEE::analyzepolarTheta -width 6
	set analyzepolarTheta 0.1
	ttk::label $w.setupTab.analyze.forceConst.polarPhi -text "phi:" -width 6 -anchor center
	ttk::entry $w.setupTab.analyze.forceConst.polarPhiEntry -textvariable ::BFEE::analyzepolarPhi -width 6
	set analyzepolarPhi 0.1
	#ttk::label $w.setupTab.analyze.forceConst.unbound -text "Unbound state:"
	#ttk::label $w.setupTab.analyze.forceConst.unrmsd -text "RMSD:" -width 6 -anchor center
	#ttk::entry $w.setupTab.analyze.forceConst.unrmsdEntry -textvariable analyzeunrmsd -width 7

	ttk::labelframe $w.setupTab.analyze.otherPar -text "Other parameter"
	ttk::label $w.setupTab.analyze.otherPar.temp -text "Temperature: "
	ttk::entry $w.setupTab.analyze.otherPar.tempEntry -textvariable ::BFEE::analyzeTemperature  -width 10
	set analyzeTemperature 300
	ttk::label $w.setupTab.analyze.otherPar.rs -text "r*: " -width 6 -anchor e
	ttk::entry $w.setupTab.analyze.otherPar.rsEntry -textvariable ::BFEE::r_star -width 10
	set r_star 30

	ttk::button $w.setupTab.analyze.compute -text "Compute Binding Free Energy" -command ::BFEE::computeFreeEnergy
    
    # drawing GUI

	grid $w.discription -columnspan 2 -sticky we
	grid $w.jobType1 $w.jobType2 -sticky we
	grid $w.setupTab -columnspan 2
	$w.setupTab add $w.setupTab.setup -text "Setup"
	$w.setupTab add $w.setupTab.analyze -text "Analyze"

	grid $w.setupTab.setup.input -sticky nswe
	grid $w.setupTab.setup.input.psfFile $w.setupTab.setup.input.psfEntry $w.setupTab.setup.input.psfButton -sticky w -padx 0.5m
	grid $w.setupTab.setup.input.coorFile $w.setupTab.setup.input.coorEntry $w.setupTab.setup.input.coorButton -sticky w -padx 0.5m
	grid $w.setupTab.setup.input.velFile $w.setupTab.setup.input.velEntry $w.setupTab.setup.input.velButton -sticky w -padx 0.5m
	grid $w.setupTab.setup.input.xscFile $w.setupTab.setup.input.xscEntry $w.setupTab.setup.input.xscButton -sticky w -padx 0.5m

	#grid $w.setupTab.setup.inputLig -sticky nswe
	#grid $w.setupTab.setup.inputLig.psfFile $w.setupTab.setup.inputLig.psfEntry $w.setupTab.setup.inputLig.psfButton -sticky w -padx 0.5m
	#grid $w.setupTab.setup.inputLig.coorFile $w.setupTab.setup.inputLig.coorEntry $w.setupTab.setup.inputLig.coorButton -sticky w -padx 0.5m
	#grid $w.setupTab.setup.inputLig.velFile $w.setupTab.setup.inputLig.velEntry $w.setupTab.setup.inputLig.velButton -sticky w -padx 0.5m
	#grid $w.setupTab.setup.inputLig.xscFile $w.setupTab.setup.inputLig.xscEntry $w.setupTab.setup.inputLig.xscButton -sticky w -padx 0.5m

	grid $w.setupTab.setup.otherPar -sticky nswe
	grid $w.setupTab.setup.otherPar.temp $w.setupTab.setup.otherPar.tempEntry -columnspan 1 -sticky w
	grid $w.setupTab.setup.otherPar.param -column 0 -rowspan 3
	grid $w.setupTab.setup.otherPar.paramDir -column 1 -row 1 -rowspan 2
	grid $w.setupTab.setup.otherPar.scroll -column 1 -row 3 -sticky nswe
	grid $w.setupTab.setup.otherPar.add -row 1 -column 2
	grid $w.setupTab.setup.otherPar.clear -row 2 -column 2
	grid $w.setupTab.setup.otherPar.sel1 $w.setupTab.setup.otherPar.sel1Entry -sticky w
	grid $w.setupTab.setup.otherPar.sel2 $w.setupTab.setup.otherPar.sel2Entry -sticky w
    grid $w.setupTab.setup.otherPar.sel3 $w.setupTab.setup.otherPar.sel3Entry -sticky w

	grid $w.setupTab.setup.generate -columnspan 2 -sticky nswe

	grid $w.setupTab.setup.contactus -sticky e
	grid $w.setupTab.setup.contactus.readme
	#grid $w.setupTab.setup.contactus.generate -columnspan 2 -sticky nswe
	#grid $w.setupTab.setup.contactus.chipot -sticky w
	#grid $w.setupTab.setup.contactus.haohao -sticky w

	grid $w.setupTab.analyze.input -sticky nswe
	grid $w.setupTab.analyze.input.boundstate -columnspan 2 -sticky w
	grid $w.setupTab.analyze.input.rmsdFile $w.setupTab.analyze.input.rmsdEntry $w.setupTab.analyze.input.rmsdButton -sticky w -padx 0.5m
	grid $w.setupTab.analyze.input.eulerThetaFile $w.setupTab.analyze.input.eulerThetaEntry $w.setupTab.analyze.input.eulerThetaButton -sticky w -padx 0.5m
	grid $w.setupTab.analyze.input.eulerPhiFile $w.setupTab.analyze.input.eulerPhiEntry $w.setupTab.analyze.input.eulerPhiButton -sticky w -padx 0.5m
	grid $w.setupTab.analyze.input.eulerPsiFile $w.setupTab.analyze.input.eulerPsiEntry $w.setupTab.analyze.input.eulerPsiButton -sticky w -padx 0.5m
	grid $w.setupTab.analyze.input.polarThetaFile $w.setupTab.analyze.input.polarThetaEntry $w.setupTab.analyze.input.polarThetaButton -sticky w -padx 0.5m
	grid $w.setupTab.analyze.input.polarPhiFile $w.setupTab.analyze.input.polarPhiEntry $w.setupTab.analyze.input.polarPhiButton -sticky w -padx 0.5m
	grid $w.setupTab.analyze.input.rFile $w.setupTab.analyze.input.rEntry $w.setupTab.analyze.input.rButton -sticky w -padx 0.5m

	grid $w.setupTab.analyze.input.unboundstate -columnspan 2 -sticky w
	grid $w.setupTab.analyze.input.unrmsdFile $w.setupTab.analyze.input.unrmsdEntry $w.setupTab.analyze.input.unrmsdButton -sticky w -padx 0.5m

	grid $w.setupTab.analyze.forceConst -sticky nswe
	grid $w.setupTab.analyze.forceConst.bound -columnspan 2 -sticky w
	grid $w.setupTab.analyze.forceConst.rmsd $w.setupTab.analyze.forceConst.rmsdEntry $w.setupTab.analyze.forceConst.eulerTheta $w.setupTab.analyze.forceConst.eulerThetaEntry $w.setupTab.analyze.forceConst.eulerPhi $w.setupTab.analyze.forceConst.eulerPhiEntry -sticky w
	grid $w.setupTab.analyze.forceConst.eulerPsi $w.setupTab.analyze.forceConst.eulerPsiEntry $w.setupTab.analyze.forceConst.polarTheta $w.setupTab.analyze.forceConst.polarThetaEntry $w.setupTab.analyze.forceConst.polarPhi $w.setupTab.analyze.forceConst.polarPhiEntry -sticky w
	#grid $w.setupTab.analyze.forceConst.unbound  -columnspan 2 -sticky w
	#grid $w.setupTab.analyze.forceConst.unrmsd $w.setupTab.analyze.forceConst.unrmsdEntry

	grid $w.setupTab.analyze.otherPar -sticky we
	grid $w.setupTab.analyze.otherPar.temp $w.setupTab.analyze.otherPar.tempEntry $w.setupTab.analyze.otherPar.rs $w.setupTab.analyze.otherPar.rsEntry -pady 1m

	grid $w.setupTab.analyze.compute -sticky we
}

# the function of Browse button
proc ::BFEE::findPsfDir {} {
	set psfTypes {
		{{Protein Structure Files} {.psf}}
		{{All Files} *}
	}
	variable psfDir
	set psfDir [tk_getOpenFile -filetypes $psfTypes]
}
proc ::BFEE::findCoorDir {} {
	set coorTypes {
		{{Coordinate Files} {.pdb .coor}}
		{{All Files} *}
	}
	variable coorDir
	set coorDir [tk_getOpenFile -filetypes $coorTypes]
}
proc ::BFEE::findVelDir {} {
	set velTypes {
		{{Velocity Files} {.vel}}
		{{All Files} *}
	}
	variable velDir
	set velDir [tk_getOpenFile -filetypes $velTypes]
}
proc ::BFEE::findXscDir {} {
	set xscTypes {
		{{eXtended System Configuration Files} {.xsc}}
		{{All Files} *}
	}
	variable xscDir
	set xscDir [tk_getOpenFile -filetypes $xscTypes]
}
proc ::BFEE::findPsfDirLig {} {
	variable psfDirLig
	set psfDirLig [tk_getOpenFile]
}
proc ::BFEE::findCoorDirLig {} {
	variable coorDirLig
	set coorDirLig [tk_getOpenFile]
}
proc ::BFEE::findVelDirLig {} {
	variable velDirLig
	set velDirLig [tk_getOpenFile]
}
proc ::BFEE::findXscDirLig {} {
	variable xscDirLig
	set xscDirLig [tk_getOpenFile]
}

# functions for add and clear button
proc ::BFEE::addParamDir {} {
	set paramTypes {
		{{Force Field Parameters Files} {.inp .prm}}
		{{All Files} *}
	}
	variable paramDirList
	set dir [tk_getOpenFile -multiple 1 -filetypes $paramTypes]
	if {$dir != ""} {
		lappend paramDirList {*}$dir
	}
}
proc ::BFEE::clearParamDir {} {
	variable paramDirList
	set paramDirList []
}

# functions for browse button in analyze notebook
set pmfTypes {
	{{Potential of Mean Force} {.pmf}}
	{{All Files} *}
}
proc ::BFEE::findrmsdDir {} {
	variable rmsdDir
	set rmsdDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findeulerThetaDir {} {
	variable eulerThetaDir
	set eulerThetaDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findeulerPhiDir {} {
	variable eulerPhiDir
	set eulerPhiDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findeulerPsiDir {} {
	variable eulerPsiDir
	set eulerPsiDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findpolarThetaDir {} {
	variable polarThetaDir
	set polarThetaDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findpolarPhiDir {} {
	variable polarPhiDir
	set polarPhiDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findrDir {} {
	variable rDir
	set rDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}
proc ::BFEE::findunrmsdDir {} {
	variable unrmsdDir
	set unrmsdDir [tk_getOpenFile -filetypes $::BFEE::pmfTypes]
}

# about us
proc ::BFEE::aboutus {} {
	variable version
	tk_messageBox -type ok -icon info \
	-message \
"Binding Free Energy Estimator $version
Simple User's Guide:
step 1: set corresponding options, generate files
step 2: run simulations, take care the steps of 
	each simulations, the <center> option of
	the harmonic restraints and the <lowerbounday>
	and the <upperboundary> of each colvar
step 3: remove the comments in each pmf file, then
	calculate binding free energy using the 
	analyze tab"
}

# move the box so that the center of mass of protein is (0,0,0)
# do it before generating input files
# generate ./eq_normalized.coor and ./eq.normalized.xst in the current folder
proc ::BFEE::move_complex {psf_file coor_file vel_file xsc_file} {
	variable sel1
	variable sel2

	mol addfile $psf_file 
	mol addfile $coor_file
	set all [atomselect top all]
	set protein [atomselect top "$sel1 and noh"]
	$all moveby [vecsub {0 0 0} [measure center $protein]]
	$all writenamdbin ./eq_normalized.coor
	pbc readxst $xsc_file
	set pbc_now [lindex [pbc get -now] 0]
	pbc set $pbc_now -now
	pbc writexst ./eq_normalized.xst
	unset all
	unset protein
	mol delete top
}

# check the polar angle theta
# do it before generating input files, if theta < 30 or > 150, then we rotate it by 90 degrees
# to avoid the rotating invarience of theta
proc ::BFEE::check_polar_theta {psf_file coor_file vel_file xsc_file} {
	variable sel1
	variable sel2

	mol addfile $psf_file 
	mol addfile $coor_file
	set selection1 [atomselect top "$sel1 and noh"]
	set selection2 [atomselect top "$sel2 and noh"]

	# calculate the angle of native state
	set c1 [measure center $selection1]
	set c2 [measure center $selection2]
	set c [vecsub $c2 $c1]
	set c [vecnorm $c]

	set angle [expr 180 / 3.1415926 * acos([lindex $c 2])]
	
	if {$angle < 30} {
		tk_messageBox -type ok -icon info -message \
"Warning: the direction vector between the protein and the ligand
is almost parallel to z-axis! Rotate automatically!" 

		# rotate the complex by 90 degrees along x axis
		# coor file
		set all [atomselect top all]
		$all move [trans center [measure center $all] axis x 90 ]
		$all writenamdbin $coor_file
		unset all
		mol delete top
		# vel file
		mol addfile $psf_file
		mol addfile $vel_file type namdbin
		set all [atomselect top all]
		$all move [trans center [measure center $all] axis x 90 ]
		$all writenamdbin $vel_file
		unset all
		#mol delete top
		# xsc file
		pbc readxst $xsc_file
		set pbc_now [lindex [pbc get -now] 0]
		# exchange y-z
		set temp [lindex $pbc_now 2]
		set pbc_now [lreplace $pbc_now 2 2 [lindex $pbc_now 1]]
		set pbc_now [lreplace $pbc_now 1 1 $temp]
		pbc set $pbc_now -now
		pbc writexst $xsc_file
		mol delete top
	} elseif {$angle > 150} {
		tk_messageBox -type ok -icon info -message \
"Warning: the direction vector between the protein and the ligand
is almost parallel to z-axis! Rotate automatically!" 

		# rotate the complex by -90 degrees along x axis
		# coor file
		set all [atomselect top all]
		$all move [trans center [measure center $all] axis x -90 ]
		$all writenamdbin $coor_file
		unset all
		mol delete top
		# vel file
		mol addfile $psf_file
		mol addfile $vel_file type namdbin
		set all [atomselect top all]
		$all move [trans center [measure center $all] axis x -90 ]
		$all writenamdbin $vel_file
		unset all
		#mol delete top
		# xsc file
		pbc readxst $xsc_file
		set pbc_now [lindex [pbc get -now] 0]
		# exchange y-z
		set temp [lindex $pbc_now 2]
		set pbc_now [lreplace $pbc_now 2 2 [lindex $pbc_now 1]]
		set pbc_now [lreplace $pbc_now 1 1 $temp]
		pbc set $pbc_now -now
		pbc writexst $xsc_file
		mol delete top
	} else {
		mol delete top
	}
}

# check the polar angle phi
# do it before generating input files, is phi > 150 or if phi < -150, then we rotate it by 180 degrees
# to avoid the periodicity of phi
proc ::BFEE::check_polar_phi {psf_file coor_file vel_file xsc_file} {
	variable sel1
	variable sel2

	mol addfile $psf_file 
	mol addfile $coor_file
	set selection1 [atomselect top "$sel1 and noh"]
	set selection2 [atomselect top "$sel2 and noh"]

	# calculate the angle of native state
	set c1 [measure center $selection1]
	set c2 [measure center $selection2]
	set c [vecsub $c2 $c1]
	set c [vecnorm $c]

	set angle [expr 180 / 3.1415926 * atan2([lindex $c 1], [lindex $c 0])]
	
	if {$angle > 150 || $angle < -150} {
		tk_messageBox -type ok -icon info -message \
"Warning: the direction vector between the protein and the ligand
is close to 180 degrees! Rotate automatically!" 

		# rotate the complex by 180 degrees along z axis
		# coor file
		set all [atomselect top all]
		$all move [trans center [measure center $all] axis z 180 ]
		$all writenamdbin $coor_file
		unset all
		mol delete top
		# vel file
		mol addfile $psf_file
		mol addfile $vel_file type namdbin
		set all [atomselect top all]
		$all move [trans center [measure center $all] axis z 180 ]
		$all writenamdbin $vel_file
		unset all
		#mol delete top
		# xsc file
		#pbc readxst $xsc_file
		#set pbc_now [lindex [pbc get -now] 0]
		# exchange y-z
		#set temp [lindex $pbc_now 2]
		#set pbc_now [lreplace $pbc_now 2 2 [lindex $pbc_now 1]]
		#set pbc_now [lreplace $pbc_now 1 1 $temp]
		#pbc set $pbc_now -now
		#pbc writexst $xsc_file
		mol delete top
	} else {
		mol delete top
	}
}

# output the definition of collective variables
proc ::BFEE::colvar_RMSD {set_boundary output_file {rmsdFile ./ligand.pdb}} {
	puts $output_file "
colvar {
    name RMSD\n"
	if {$set_boundary == 1} {
		puts $output_file "
    width 0.05 
    lowerboundary 0 
    upperboundary 3

    subtractAppliedForce on
    expandboundaries  on
    extendedLagrangian on
    extendedFluctuation 0.05\n"
	}

	puts $output_file "
    rmsd {
        atoms {
            atomsFile $rmsdFile
	    atomsCol B
	    atomsColValue 1.0
        }
        refpositionsfile  $rmsdFile
    }
}\n"
}

proc ::BFEE::colvar_euler_angle {set_boundary output_file euler_angle} {
	puts $output_file "
colvar {
    name $euler_angle
    scriptedFunction  $euler_angle \n"
	if {$set_boundary == 1} {
		puts $output_file "
    width                 1 
                       
    lowerboundary         -10 
    upperboundary          10

    subtractAppliedForce on
    expandboundaries  on
    extendedLagrangian      on
    extendedFluctuation    1\n"
	}

	puts $output_file "
    Orientation {
        atoms { 
            atomsFile ./ligand.pdb
	    atomsCol B
	    atomsColValue 1.0
            centerReference    on
            rotateReference   on
	    enableFitGradients no
            fittingGroup {
                atomsFile ./protein.pdb
	        atomsCol B
	        atomsColValue 1.0
            }   
            refPositionsFile  ./protein.pdb
         }   
         refPositionsFile  ./ligand.pdb
    }   
}\n"
}

# calculate the value of polar angles
# will be used in polar angle CVs
proc ::BFEE::polar_angle_value {polar_angle psf_file coor_file} {
    variable sel1
	variable sel2

	mol addfile $psf_file 
	mol addfile $coor_file
	set selection1 [atomselect top "$sel1 and noh"]
	set selection2 [atomselect top "$sel2 and noh"]

	# calculate the angle of native state
	set c1 [measure center $selection1]
	set c2 [measure center $selection2]
	set c [vecsub $c2 $c1]
	set c [vecnorm $c]

	if {$polar_angle == "polarTheta"} {
		set angle_center [expr 180 / 3.1415926 * acos([lindex $c 2])]
	} elseif {$polar_angle == "polarPhi"} {
		set angle_center [expr 180 / 3.1415926 * atan2([lindex $c 1], [lindex $c 0])]
	}
    mol delete top
	return [expr floor($angle_center)]
}

proc ::BFEE::colvar_polar_angle {set_boundary output_file polar_angle value} {

	puts $output_file "
colvar {
    name $polar_angle\n"
	if {$set_boundary == 1} {
		puts $output_file "
    width                 1 
                       
    lowerboundary         [expr $value - 10]
    upperboundary         [expr $value + 10]

    subtractAppliedForce    on
    expandboundaries  on
    extendedLagrangian      on
    extendedFluctuation      1\n"
	}

	puts $output_file "
    $polar_angle {
        atoms { 
            atomsFile ./ligand.pdb
	    atomsCol B
	    atomsColValue 1.0
        centerReference   on
        rotateReference   on
        fittingGroup {
            atomsFile ./protein.pdb
	    atomsCol B
	    atomsColValue 1.0
        }   
        refPositionsFile  ./protein.pdb
        }
    }
}\n"
}

# calculate the value of r*
# will be used in the distance CV 
proc ::BFEE::r_value {psf_file coor_file} {
    variable sel1
	variable sel2

	mol addfile $psf_file 
	mol addfile $coor_file
	set selection1 [atomselect top "$sel1 and noh"]
	set selection2 [atomselect top "$sel2 and noh"]

	# calculate r of native state
	set c1 [measure center $selection1]
	set c2 [measure center $selection2]

	set distance [vecdist $c1 $c2]
    return $distance
}

proc ::BFEE::colvar_r {set_boundary output_file distance} {

	set upper_distance [expr floor($distance + 21)]
	set lower_distance [expr floor($distance - 2)]

	if {$lower_distance < 0} {
		set lower_distance 0.2
	}

	puts $output_file "
colvar {
    name r\n"
    
    if {$set_boundary == 1} {
		puts $output_file "
    width                0.1 
                       
    lowerboundary        $lower_distance
    upperboundary        $upper_distance
                       
    subtractAppliedForce on
    expandboundaries  on
    extendedLagrangian      on
    extendedFluctuation    0.1\n"
    }
    
    puts $output_file "
    distance {
        forceNoPBC       yes
        group1 { 
            atomsFile ./protein.pdb
	    atomsCol B
	    atomsColValue 1.0
	}
        group2 {
            atomsFile ./ligand.pdb
	    atomsCol B
	    atomsColValue 1.0
        }
    }
}\n"
}

# output the head of colvars config file
proc ::BFEE::colvar_head {output_file} {
	puts $output_file "colvarsTrajFrequency      500
colvarsRestartFrequency   50000
"
}

# output the bias
proc ::BFEE::bias_walls {output_file colvar lowerwall upperwall} {
    puts $output_file "
harmonicWalls {
  colvars $colvar
  lowerWalls $lowerwall
  upperWalls $upperwall
  lowerWallConstant 0.2
  upperWallConstant 0.2
}\n"
}

# ti = 0 for plain harmonic restrants, 1 for forward Ti calculation, 2 for backward Ti calculation
proc ::BFEE::bias_harmonic {output_file colvar constant center {ti 0}} {
	puts $output_file "
harmonic {
   colvars      $colvar
   forceConstant   $constant
   centers  $center\n"
   if {$ti != 0} {
       puts $output_file "targetNumSteps 2000000
   targetEquilSteps 100000
   targetForceConstant 0
   targetForceExponent 4\n"
   }
   if {$ti == 1} {
       puts $output_file "lambdaSchedule 1 0.99999 0.9999 0.999 0.99 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0\n"
   }
   if {$ti == 2} {
       puts $output_file "lambdaSchedule 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.999 0.9999 0.99999 1\n"
   }
   puts $output_file "}\n"
}

proc ::BFEE::bias_abf {output_file colvar} {
	puts $output_file "
abf {
   colvars           $colvar
   FullSamples       10000
   historyfreq       50000
   writeCZARwindowFile
}
metadynamics {
   colvars           $colvar
   hillWidth         3.0
   hillWeight        0.05
   wellTempered      on
   biasTemperature   4000
}\n"
}

# output the colvar that restrains the protein
proc ::BFEE::colvar_protein {output_file psf_file coor_file ref_protein} {
	variable sel1

	mol addfile $psf_file 
	mol addfile $coor_file
	set selection1 [atomselect top "$sel1 and noh"]

	puts $output_file "
colvar { 
  name translation 
  distance { 
    group1 { 
      atomsFile ./protein.pdb
      atomsCol B
      atomsColValue 1.0 
    } 
    group2 { 
      dummyAtom ([lindex [measure center $selection1] 0], [lindex [measure center $selection1] 1], [lindex [measure center $selection1] 2])
    } 
  } 
} 
harmonic { 
  colvars translation 
  centers 0.0 
  forceConstant 10.0 
} 


colvar { 
  name orientation 
  orientation { 
    atoms { 
      atomsFile ./protein.pdb
      atomsCol B
      atomsColValue 1.0 
    } 
    refPositionsFile   $ref_protein
  } 
} 
harmonic { 
  colvars orientation 
  centers (1.0, 0.0, 0.0, 0.0)
  forceConstant 200.0 
} \n"

	mol delete top
}

# write reference files for colvars
proc ::BFEE::write_ref {psf_file coor_file} {
	variable sel1
	variable sel2
	
	mol addfile $psf_file 
	mol addfile $coor_file

	set all [atomselect top all]
	set selection1 [atomselect top "$sel1 and noh"]
	set selection2 [atomselect top "$sel2 and noh"]
	$all set beta 0
	$selection1 set beta 1
	$all writepdb protein.pdb
	$all set beta 0
	$selection2 set beta 1
	$all writepdb ligand.pdb

	$all delete
	mol delete top
}

# add water box based on psf and coor (for step 7)
proc ::BFEE::add_watbox {psf_file coor_file} {
	variable sel1
	variable sel2

	mol addfile $psf_file
	mol addfile $coor_file

	set selection1 [atomselect top "$sel1 and noh"]
	set selection2 [atomselect top "$sel2 and noh"]

	set direction [vecsub [measure center $selection2] [measure center $selection1]]

	# get the direction of protein-ligand distance vector
	# for extending water box
	if {[lindex $direction 0] > 0} {
		set xm 0
		set xp 22
	} else {
		set xp 0
		set xm 22
	}
	if {[lindex $direction 1] > 0} {
		set ym 0
		set yp 22
	} else {
		set yp 0
		set ym 22
	}
	if {[lindex $direction 2] > 0} {
		set zm 0
		set zp 22
	} else {
		set zp 0
		set zm 22
	}

	set all [atomselect top all]
	$all writepdb $coor_file.pdb
	mol delete top

	#solvate $psf_file $coor_file.pdb -t 15 -o solvated -s W2 -b 2.2
	solvate $psf_file $coor_file.pdb -x $xm -y $ym -z $zm +x $xp +y $yp +z $zp -o solvated -s W2 -b 2.2
	mol delete top
}

# remove protein based on psf and coor (for step 8)
proc ::BFEE::remove_protein {psf_file coor_file} {
	variable sel3
	variable sel2

	mol addfile $psf_file
	mol addfile $coor_file
	set nopro [atomselect top "not $sel3"]
	$nopro writepsf lig.psf
	$nopro writepdb lig.pdb

	$nopro delete
	mol delete top
}

# output PBC information for namd config file (for step 7,8)
proc ::BFEE::write_cell {output_file temperature pdb_file} {

	mol new $pdb_file

	puts $output_file "coordinates $pdb_file"
	puts $output_file "temperature $temperature"

	set all [atomselect top all] 
	set minmax [measure minmax $all] 
	set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]] 
	puts $output_file "cellBasisVector1 [lindex $vec 0] 0 0" 
	puts $output_file "cellBasisVector2 0 [lindex $vec 1] 0" 
	puts $output_file "cellBasisVector3 0 0 [lindex $vec 2]" 
	set center [measure center $all] 
	puts $output_file "cellOrigin $center" 
	
	mol delete top
}

# output namd config file
proc ::BFEE::namd_config {output_file psf_file restart restart_prefix namd_output param temperature numsteps pdb_file {cv_file colvar.in} {fep 0} {fepfile ./fep.fep}} {
	puts $output_file "
##################################################
# PRODUCTION MD
##################################################


##################################################
# MD SECTION
##################################################



# TOPOLOGY

structure            $psf_file


# FORCE FIELD
"
	foreach i $param {
	puts $output_file "parameters       ../[file tail $i]"
}

	puts $output_file "
paraTypeCharmm       on


# 1-4 TERMs

exclude              scaled1-4
1-4scaling           1.0


# INPUT FILES
"
	if {$restart == 1} {
		puts $output_file "
coordinates          $pdb_file
bincoordinates       $restart_prefix.coor
binvelocities        $restart_prefix.vel
ExtendedSystem       $restart_prefix.xsc \n"
	} else {
		write_cell $output_file $temperature $restart_prefix
	}

	puts $output_file "
# OUTPUT FILES

binaryoutput         yes
binaryrestart        yes

outputname           output/$namd_output
restartname          output/$namd_output


# DCD FILE

dcdFile              output/$namd_output.dcd
dcdUnitCell          yes


# FREQUENCY FOR DUMPING OUTPUT DATA

outputenergies       5000
outputtiming         5000
outputpressure       5000
restartfreq          5000
XSTFreq              5000
dcdFreq              5000


# CUT-OFFs

hgroupcutoff         2.8
switching            on
switchdist           10.0
cutoff               12.0
pairlistdist         14.0


# WRAPPING

wrapAll              off
wrapWater            on


# CONSTANT-T

stochRescale             on
stochRescaleTemp         $temperature
stochRescalePeriod       2.0 


# CONSTANT-P

langevinpiston       on 
langevinpistontarget 1.01325 
langevinpistonperiod 200
langevinpistondecay  100
langevinpistontemp   $temperature

strainrate           0.  0.  0.
usegrouppressure     yes

useflexiblecell      no  


# PME

PME                  yes
PMETolerance         10e-6
PMEInterpOrder       4
PMEGridSpacing       1.0


# MULTIPLE TIME-STEP PROPAGATOR

timestep             2.0

fullelectfrequency   2
nonbondedfreq        1


# SHAKE/RATTLE

rigidbonds           all     
rigidtolerance       0.00001
rigiditerations      400


# PARALLELISM

stepspercycle        10
splitpatch           hydrogen
margin               2


##################################################
# ABF SECTION
##################################################

colvars              on
colvarsConfig        $cv_file

source ../CVs.tcl\n"

    if {$fep == 0} {
        puts $output_file "
minimize             100
reinitvels           $temperature

run $numsteps\n"
    } else {
        puts $output_file "
source fep.tcl
alch on
alchType FEP
alchFile $fepfile
alchCol B
alchOutFile $namd_output.fepout
alchOutFreq 10
alchVdwLambdaEnd 1.0
alchElecLambdaStart 0.5
alchEquilSteps 100000\n"
    }
    if {$fep == 1} {
        puts $output_file "
runFEPmin 0.0 1.0 0.02 500000 1000 $temperature\n
"
    }
    if {$fep == 2} {
        puts $output_file "
runFEPmin 1.0 0.0 -0.02 500000 1000 $temperature\n
"
    }
}

# output namd config file for the first restart
# fep = 0 for equilibration, 1 for creation, 2 for annihilation
proc ::BFEE::namd_restart {output_file psf_file restart_prefix namd_output param temperature numsteps pdb_file {cv_file colvar.in}} {
	puts $output_file "
##################################################
# PRODUCTION MD
##################################################


##################################################
# MD SECTION
##################################################



# TOPOLOGY

structure            $psf_file


# FORCE FIELD
"
	foreach i $param {
	puts $output_file "parameters       ../[file tail $i]"
}

	puts $output_file "
paraTypeCharmm       on


# 1-4 TERMs

exclude              scaled1-4
1-4scaling           1.0


# INPUT FILES
coordinates          $pdb_file
bincoordinates       $restart_prefix.coor
binvelocities        $restart_prefix.vel
ExtendedSystem       $restart_prefix.xsc


# OUTPUT FILES

binaryoutput         yes
binaryrestart        yes

outputname           output/$namd_output
restartname          output/$namd_output


# DCD FILE

dcdFile              output/$namd_output.dcd
dcdUnitCell          yes


# FREQUENCY FOR DUMPING OUTPUT DATA

outputenergies       5000
outputtiming         5000
outputpressure       5000
restartfreq          5000
XSTFreq              5000
dcdFreq              5000


# CUT-OFFs

hgroupcutoff         2.8
switching            on
switchdist           10.0
cutoff               12.0
pairlistdist         14.0


# WRAPPING

wrapAll              off
wrapWater            on


# CONSTANT-T

stochRescale             on
stochRescaleTemp         $temperature
stochRescalePeriod       2.0 


# CONSTANT-P

langevinpiston       on 
langevinpistontarget 1.01325
langevinpistonperiod 200
langevinpistondecay  100
langevinpistontemp   $temperature

strainrate           0.  0.  0.
usegrouppressure     yes

useflexiblecell      no  


# PME

PME                  yes
PMETolerance         10e-6
PMEInterpOrder       4
PMEGridSpacing       1.0


# MULTIPLE TIME-STEP PROPAGATOR

timestep             2.0

fullelectfrequency   2
nonbondedfreq        1


# SHAKE/RATTLE

rigidbonds           all     
rigidtolerance       0.00001
rigiditerations      400


# PARALLELISM

stepspercycle        10
splitpatch           hydrogen
margin               2


##################################################
# ABF SECTION
##################################################

colvars              on
colvarsConfig        colvar.in
colvarsInput         $restart_prefix.colvars.state

source ../CVs.tcl

run $numsteps"
}

# output CVs.tcl
proc ::BFEE::writeCV {output_file} {
	puts $output_file "# Euler angles
# Phi
namespace eval eulerPhi { } 
proc calc_eulerPhi { args } {

    global eulerPhi::q0
    global eulerPhi::q1
    global eulerPhi::q2
    global eulerPhi::q3

    set q0 \[ lindex \[ lindex \$args 0 \] 0 \]
    set q1 \[ lindex \[ lindex \$args 0 \] 1 \]
    set q2 \[ lindex \[ lindex \$args 0 \] 2 \]
    set q3 \[ lindex \[ lindex \$args 0 \] 3 \]

    set f  \[ expr 180 / 3.1415926 * atan2(2 * (\$q0 * \$q1 + \$q2 * \$q3), 1 - 2 * (\$q1 * \$q1 + \$q2 * \$q2)) \]

    return \$f
}


proc calc_eulerPhi_gradient { args } {

    global eulerPhi::q0
    global eulerPhi::q1
    global eulerPhi::q2
    global eulerPhi::q3

    set gq0 \[expr 180 / 3.1415926 * 2 * \$q1 * (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1) / ((2 * \$q0 * \$q1 + 2 * \$q2 * \$q3)**2 + (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1)**2)\]
    set gq1 \[expr 180 / 3.1415926 * (2 * \$q0 * (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1) - 4 * \$q1 * (-2 * \$q0 * \$q1 - 2 * \$q2 * \$q3)) / ((2 * \$q0 * \$q1 + 2 * \$q2 * \$q3)**2 + (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1)**2)\]
    set gq2 \[expr 180 / 3.1415926 * (-4 * \$q2 * (-2 * \$q0 * \$q1 - 2 * \$q2 * \$q3) + 2 * \$q3 * (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1)) / ((2 * \$q0 * \$q1 + 2 * \$q2 * \$q3)**2 + (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1)**2)\]
    set gq3 \[expr 180 / 3.1415926 * 2 * \$q2 * (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1) / ((2 * \$q0 * \$q1 + 2 * \$q2 * \$q3)**2 + (-2 * \$q1 * \$q1 - 2 * \$q2 * \$q2 + 1)**2)\]

    set g \"{ \$gq0 \$gq1 \$gq2 \$gq3 }\"

    return \$g
}


# Psi
namespace eval eulerPsi { } 
proc calc_eulerPsi { args } {

    global eulerPsi::q0
    global eulerPsi::q1
    global eulerPsi::q2
    global eulerPsi::q3

    set q0 \[ lindex \[ lindex \$args 0 \] 0 \]
    set q1 \[ lindex \[ lindex \$args 0 \] 1 \]
    set q2 \[ lindex \[ lindex \$args 0 \] 2 \]
    set q3 \[ lindex \[ lindex \$args 0 \] 3 \]

    set f  \[ expr 180 / 3.1415926 * atan2(2 * (\$q0 * \$q3 + \$q1 * \$q2), 1 - 2 * (\$q2 * \$q2 + \$q3 * \$q3)) \]

    return \$f
}

proc calc_eulerPsi_gradient { args } {

    global eulerPsi::q0
    global eulerPsi::q1
    global eulerPsi::q2
    global eulerPsi::q3

    set gq0 \[expr 180 / 3.1415926 * 2 * \$q3 * (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1) / ((2 * \$q0 * \$q3 + 2 * \$q1 * \$q2)**2 + (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1)**2)\]
    set gq1 \[expr 180 / 3.1415926 * 2 * \$q2 * (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1) / ((2 * \$q0 * \$q3 + 2 * \$q1 * \$q2)**2 + (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1)**2)\]
    set gq2 \[expr 180 / 3.1415926 * (2 * \$q1 * (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1) - 4 * \$q2 * (-2 * \$q0 * \$q3 - 2 * \$q1 * \$q2)) / ((2 * \$q0 * \$q3 + 2 * \$q1 * \$q2)**2 + (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1)**2)\]
    set gq3 \[expr 180 / 3.1415926 * (2 * \$q0 * (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1) - 4 * \$q3 * (-2 * \$q0 * \$q3 - 2 * \$q1 * \$q2)) / ((2 * \$q0 * \$q3 + 2 * \$q1 * \$q2)**2 + (-2 * \$q2 * \$q2 - 2 * \$q3 * \$q3 + 1)**2)\]

    set g \"{ \$gq0 \$gq1 \$gq2 \$gq3 }\"

    return \$g
}

# Theta
namespace eval eulerTheta { } 
proc calc_eulerTheta { args } {

    global eulerTheta::q0
    global eulerTheta::q1
    global eulerTheta::q2
    global eulerTheta::q3

    set q0 \[ lindex \[ lindex \$args 0 \] 0 \]
    set q1 \[ lindex \[ lindex \$args 0 \] 1 \]
    set q2 \[ lindex \[ lindex \$args 0 \] 2 \]
    set q3 \[ lindex \[ lindex \$args 0 \] 3 \]

    set f  \[ expr 180 / 3.1415926 * asin(2 * (\$q0 * \$q2 - \$q3 * \$q1)) \]

    return \$f
}

proc calc_eulerTheta_gradient { args } {

    global eulerTheta::q0
    global eulerTheta::q1
    global eulerTheta::q2
    global eulerTheta::q3

    set gq0 \[expr 180 / 3.1415926 * 2 * \$q2 / sqrt(-((2 * \$q0 * \$q2 - 2 * \$q1 * \$q3)**2) + 1)\]
    set gq1 \[expr 180 / 3.1415926 * -2 * \$q3 / sqrt(-((2 * \$q0 * \$q2 - 2 * \$q1 * \$q3)**2) + 1)\]
    set gq2 \[expr 180 / 3.1415926 * 2 * \$q0 / sqrt(-((2 * \$q0 * \$q2 - 2 * \$q1 * \$q3)**2) + 1)\]
    set gq3 \[expr 180 / 3.1415926 * -2 * \$q1 / sqrt(-((2 * \$q0 * \$q2 - 2 * \$q1 * \$q3)**2) + 1)\]

    set g \"{ \$gq0 \$gq1 \$gq2 \$gq3 }\"

    return \$g
}"
}

# generate directory
proc ::BFEE::gen_subdir {name} {
	file mkdir $name
	cd $name 
	file mkdir output
	cd ..
}
proc ::BFEE::gen_dir {} {
	gen_subdir 001_RMSD_bound
	gen_subdir 002_euler_theta
	gen_subdir 003_euler_phi
	gen_subdir 004_euler_psi
	gen_subdir 005_polar_theta
	gen_subdir 006_polar_phi
	gen_subdir 007_r
	gen_subdir 008_RMSD_unbound
    gen_subdir alter_alchemy
}

# copy files to the working directory
proc ::BFEE::copy_files {} {
	variable psfDir
	variable coorDir
	variable velDir
	variable xscDir
	variable psfDirLig
	variable coorDirLig
	variable velDirLig
	variable xscDirLig
	variable paramDirList
	variable sel1
	variable sel2
    variable sel3

	check_polar_theta $psfDir $coorDir $velDir $xscDir
	check_polar_phi $psfDir $coorDir $velDir $xscDir
	move_complex $psfDir $coorDir $velDir $xscDir

	file copy -force $psfDir ./bound.psf
	file copy -force ./eq_normalized.coor ./eq.coor
	file copy -force $velDir ./eq.vel
	file copy -force ./eq_normalized.xst ./eq.xsc

	file copy -force ./bound.psf ./007_r/bound.psf
	file copy -force ./eq_normalized.coor ./007_r/bound.coor

	file copy -force ./bound.psf ./008_RMSD_unbound/bound.psf
	file copy -force ./eq_normalized.coor ./008_RMSD_unbound/bound.coor

	foreach i $paramDirList {
		file copy -force $i ./
	}

	mol addfile $psfDir 
	mol addfile ./eq_normalized.coor
	set all [atomselect top all]
	$all writepdb ./bound.pdb
    
    # write fep file
    cd alter_alchemy
    $all set beta 0
    set lig [atomselect top $sel2]
    $lig set beta 1
    $all writepdb fep.fep
    set all_Except_Lig [atomselect top "not $sel3"]
    $all_Except_Lig writepsf unbound.psf
    $all_Except_Lig writepdb unbound.pdb
    $all_Except_Lig writepdb fep_unbound.fep
    $all set beta 0
    set lig_noh [atomselect top "$sel2 and noh"]
    $lig_noh set beta 1
    $all_Except_Lig writepdb ligand_unbound.pdb
    cd ..
    
	mol delete top
}

# generate colvar.in in each directory
proc ::BFEE::write_colvar_files {} {

    # calculate the value of polarTheta, polarPhi and R
    set pTheta [polar_angle_value polarTheta ./bound.psf ./bound.pdb]
	set pPhi [polar_angle_value polarPhi ./bound.psf ./bound.pdb]
    set rR [r_value ./bound.psf ./bound.pdb]
    set urR [expr floor($rR + 21)]
	set lrR [expr floor($rR - 2)]

	if {$lrR < 0} {
		set lrR 0.2
	}
    
	cd 001_RMSD_bound
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 1 $of
    bias_walls $of RMSD 0 3
	bias_abf $of RMSD
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
	cd ..

	cd 002_euler_theta
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 1 $of eulerTheta
    #bias_walls $of RMSD 0 3
    bias_walls $of eulerTheta -10 10
	bias_abf $of eulerTheta
	bias_harmonic $of RMSD 10 0
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
	cd ..

	cd 003_euler_phi
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 1 $of eulerPhi
    #bias_walls $of RMSD 0 3
    #bias_walls $of eulerTheta -10 10
    bias_walls $of eulerPhi -10 10
	bias_abf $of eulerPhi
	bias_harmonic $of eulerTheta 0.1 0
	bias_harmonic $of RMSD 10 0
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
	cd ..

	cd 004_euler_psi
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 1 $of eulerPsi
    #bias_walls $of RMSD 0 3
    #bias_walls $of eulerTheta -10 10
    #bias_walls $of eulerPhi -10 10
    bias_walls $of eulerPsi -10 10
	bias_abf $of eulerPsi
	bias_harmonic $of eulerPhi 0.1 0 
	bias_harmonic $of eulerTheta 0.1 0
	bias_harmonic $of RMSD 10 0
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
	cd ..

	cd 005_polar_theta
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 0 $of eulerPsi
	colvar_polar_angle 1 $of polarTheta $pTheta
    #bias_walls $of RMSD 0 3
    #bias_walls $of eulerTheta -10 10
    #bias_walls $of eulerPhi -10 10
    #bias_walls $of eulerPsi -10 10
    bias_walls $of polarTheta [expr $pTheta - 10] [expr $pTheta + 10]
	bias_abf $of polarTheta
	bias_harmonic $of eulerPsi 0.1 0
	bias_harmonic $of eulerPhi 0.1 0 
	bias_harmonic $of eulerTheta 0.1 0
	bias_harmonic $of RMSD 10 0
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
	cd ..

	cd 006_polar_phi
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 0 $of eulerPsi
    colvar_polar_angle 0 $of polarTheta $pTheta
	colvar_polar_angle 1 $of polarPhi $pPhi
    #bias_walls $of RMSD 0 3
    #bias_walls $of eulerTheta -10 10
    #bias_walls $of eulerPhi -10 10
    #bias_walls $of eulerPsi -10 10
    #bias_walls $of polarTheta [expr $pTheta - 10] [expr $pTheta + 10]
    bias_walls $of polarPhi [expr $pPhi - 10] [expr $pPhi + 10]
	bias_abf $of polarPhi
	bias_harmonic $of polarTheta 0.1 $pTheta
	bias_harmonic $of eulerPsi 0.1 0
	bias_harmonic $of eulerPhi 0.1 0 
	bias_harmonic $of eulerTheta 0.1 0
	bias_harmonic $of RMSD 10 0
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
	cd ..

	cd 007_r
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 0 $of eulerPsi
    colvar_polar_angle 0 $of polarTheta $pTheta
	colvar_polar_angle 0 $of polarPhi $pPhi
	colvar_r 1 $of $rR
    #bias_walls $of RMSD 0 3
    #bias_walls $of eulerTheta -10 10
    #bias_walls $of eulerPhi -10 10
    #bias_walls $of eulerPsi -10 10
    #bias_walls $of polarTheta [expr $pTheta - 10] [expr $pTheta + 10]
    #bias_walls $of polarPhi [expr $pPhi - 10] [expr $pPhi + 10]
    bias_walls $of r $lrR $urR
	bias_abf $of r
	bias_harmonic $of polarPhi 0.1 $pPhi
	bias_harmonic $of polarTheta 0.1 $pTheta
	bias_harmonic $of eulerPsi 0.1 0
	bias_harmonic $of eulerPhi 0.1 0 
	bias_harmonic $of eulerTheta 0.1 0
	bias_harmonic $of RMSD 10 0
	colvar_protein $of ./solvated.psf ./solvated.pdb ./solvated.pdb
	close $of
	cd ..

	cd 008_RMSD_unbound
	set of [open colvar.in w]
	colvar_head $of
	colvar_RMSD 1 $of
    bias_walls $of RMSD 0 3
	bias_abf $of RMSD
	close $of
	cd ..
    
    cd alter_alchemy
    # 4 .in files
    set of [open unbound_genAtom.in w]
    colvar_head $of
	colvar_RMSD 0 $of ./ligand_unbound.pdb
    bias_harmonic $of RMSD 10 0
	close $of
    
    set of [open unbound_genRestraint.in w]
    colvar_head $of
	colvar_RMSD 0 $of ./ligand_unbound.pdb
    bias_harmonic $of RMSD 10 0 1
	close $of
    
    set of [open unbound_annihRestraint.in w]
    colvar_head $of
	colvar_RMSD 0 $of ./ligand_unbound.pdb
    bias_harmonic $of RMSD 10 0 2
	close $of
    
    set of [open bound_genAtom.in w]
    colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 0 $of eulerPsi
    colvar_polar_angle 0 $of polarTheta $pTheta
	colvar_polar_angle 0 $of polarPhi $pPhi
	colvar_r 0 $of $rR
    bias_harmonic $of RMSD 10 0
	bias_harmonic $of r 10 $rR
	bias_harmonic $of polarPhi 0.1 $pPhi
	bias_harmonic $of polarTheta 0.1 $pTheta
	bias_harmonic $of eulerPsi 0.1 0
	bias_harmonic $of eulerPhi 0.1 0 
	bias_harmonic $of eulerTheta 0.1 0
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
    
    set of [open bound_genRestraint.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 0 $of eulerPsi
    colvar_polar_angle 0 $of polarTheta $pTheta
	colvar_polar_angle 0 $of polarPhi $pPhi
	colvar_r 0 $of $rR
    bias_harmonic $of RMSD 10 0 1
	bias_harmonic $of r 10 $rR 1
	bias_harmonic $of polarPhi 0.1 $pPhi 1
	bias_harmonic $of polarTheta 0.1 $pTheta 1
	bias_harmonic $of eulerPsi 0.1 0 1
	bias_harmonic $of eulerPhi 0.1 0  1
	bias_harmonic $of eulerTheta 0.1 0 1
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
    
    set of [open bound_annihRestraint.in w]
	colvar_head $of
	colvar_RMSD 0 $of
	colvar_euler_angle 0 $of eulerTheta
	colvar_euler_angle 0 $of eulerPhi
	colvar_euler_angle 0 $of eulerPsi
    colvar_polar_angle 0 $of polarTheta $pTheta
	colvar_polar_angle 0 $of polarPhi $pPhi
	colvar_r 0 $of $rR
    bias_harmonic $of RMSD 10 0 2
	bias_harmonic $of r 10 $rR 2
	bias_harmonic $of polarPhi 0.1 $pPhi 2
	bias_harmonic $of polarTheta 0.1 $pTheta 2
	bias_harmonic $of eulerPsi 0.1 0 2
	bias_harmonic $of eulerPhi 0.1 0  2
	bias_harmonic $of eulerTheta 0.1 0 2
	colvar_protein $of ../bound.psf ../eq.coor ./protein.pdb
	close $of
    
    cd ..
}

# generate namd config files in each directory
proc ::BFEE::write_namd_files {} {
	variable paramDirList
	variable temperature

	set cv [open CVs.tcl w]
	writeCV $cv
	close $cv

	cd 001_RMSD_bound
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ../bound.psf 1 ../eq RMSD_bound $paramDirList $temperature 20000000 ../bound.pdb
	namd_restart $of2 ../bound.psf output/RMSD_bound RMSD_bound2 $paramDirList $temperature 20000000 ../bound.pdb
	close $of
	close $of2
	write_ref ../bound.psf ../eq.coor
	cd ..

	cd 002_euler_theta
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ../bound.psf 1 ../eq euler_theta $paramDirList $temperature 4000000 ../bound.pdb
	namd_restart $of2 ../bound.psf output/euler_theta euler_theta2 $paramDirList $temperature 4000000 ../bound.pdb
	close $of
	close $of2
	write_ref ../bound.psf ../eq.coor
	cd ..

	cd 003_euler_phi
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ../bound.psf 1 ../eq euler_phi $paramDirList $temperature 4000000 ../bound.pdb
	namd_restart $of2 ../bound.psf output/euler_phi euler_phi2 $paramDirList $temperature 4000000 ../bound.pdb
	close $of
	close $of2
	write_ref ../bound.psf ../eq.coor
	cd ..

	cd 004_euler_psi
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ../bound.psf 1 ../eq euler_psi $paramDirList $temperature 4000000 ../bound.pdb
	namd_restart $of2 ../bound.psf output/euler_psi euler_psi2 $paramDirList $temperature 4000000 ../bound.pdb
	close $of
	close $of2
	write_ref ../bound.psf ../eq.coor
	cd ..

	cd 005_polar_theta
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ../bound.psf 1 ../eq polar_theta $paramDirList $temperature 4000000 ../bound.pdb
	namd_restart $of2 ../bound.psf output/polar_theta polar_theta2 $paramDirList $temperature 4000000 ../bound.pdb
	close $of
	close $of2
	write_ref ../bound.psf ../eq.coor
	cd ..

	cd 006_polar_phi
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ../bound.psf 1 ../eq polar_phi $paramDirList $temperature 4000000 ../bound.pdb
	namd_restart $of2 ../bound.psf output/polar_phi polar_phi2 $paramDirList $temperature 4000000 ../bound.pdb
	close $of
	close $of2
	write_ref ../bound.psf ../eq.coor
	cd ..

	cd 007_r
	add_watbox ./bound.psf ./bound.coor
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ./solvated.psf 0 ./solvated.pdb r $paramDirList $temperature 80000000 nopdb
	namd_restart $of2 ./solvated.psf output/r r2 $paramDirList $temperature 80000000 ./solvated.pdb
	close $of
	close $of2
	write_ref ./solvated.psf ./solvated.pdb
	cd ..

	cd 008_RMSD_unbound
	remove_protein ./bound.psf ./bound.coor
	set of [open abf.conf w]
	set of2 [open abf_restart.conf w]
	namd_config $of ./lig.psf 0 ./lig.pdb RMSD_unbound $paramDirList $temperature 80000000 nopdb
	namd_restart $of2 ./lig.psf output/RMSD_unbound RMSD_unbound2 $paramDirList $temperature 80000000 ./lig.pdb
	close $of
	close $of2
	write_ref ./lig.psf ./lig.pdb
	cd ..
    
    cd alter_alchemy
    set of [open unbound_genAtom.conf w]
    namd_config $of ./unbound.psf 0 ./unbound.pdb unbound_genAtom $paramDirList $temperature 0 ./unbound.pdb ./unbound_genAtom.in 1 ./fep_unbound.fep
    close $of
    set of [open bound_genAtom.conf w]
    namd_config $of ../bound.psf 0 ../bound.pdb bound_genAtom $paramDirList $temperature 0 ../bound.pdb ./bound_genAtom.in 1 ./fep.fep
    close $of
    set of [open unbound_annihAtom.conf w]
    namd_config $of ./unbound.psf 0 ./unbound.pdb unbound_annihAtom $paramDirList $temperature 0 ./unbound.pdb ./unbound_genAtom.in 2 ./fep_unbound.fep
    close $of
    set of [open bound_annihAtom.conf w]
    namd_config $of ../bound.psf 0 ../bound.pdb bound_annihAtom $paramDirList $temperature 0 ../bound.pdb ./bound_genAtom.in 2 ./fep.fep
    close $of
    set of [open unbound_genRestraint.conf w]
    namd_config $of ./unbound.psf 0 ./unbound.pdb unbound_genRestraint $paramDirList $temperature 30000000 ./unbound.pdb ./unbound_genRestraint.in
    close $of
    set of [open bound_genRestraint.conf w]
    namd_config $of ../bound.psf 0 ../bound.pdb bound_genRestraint $paramDirList $temperature 30000000 ../bound.pdb ./unbound_genRestraint.in
    close $of
    set of [open unbound_annihRestraint.conf w]
    namd_config $of ./unbound.psf 0 ./unbound.pdb unbound_annihRestraint $paramDirList $temperature 30000000 ./unbound.pdb ./unbound_annihRestraint.in
    close $of
    set of [open bound_annihRestraint.conf w]
    namd_config $of ../bound.psf 0 ../bound.pdb bound_annihRestraint $paramDirList $temperature 30000000 ../bound.pdb ./unbound_annihRestraint.in
    close $of
    write_ref ../bound.psf ../eq.coor
    cd ..
}

# generate input files
proc ::BFEE::generateFiles {} {
	gen_dir
	copy_files
	write_namd_files
	write_colvar_files
	tk_messageBox -type ok -icon info \
	-message \
"Input files have been successfully generated!" 
}

proc ::BFEE::beta {} {
	variable analyzeTemperature
	variable BOLTZMANN
	return [expr 1.0 / ($BOLTZMANN * $analyzeTemperature)]
}

proc ::BFEE::splitline {line} {
	set line2 [string trim $line]
	regsub -all {[[:blank:]]+} $line2 " " line3
	return [split $line3]
}

# read namd pmf file
proc ::BFEE::readPMF {pmfFile} {

	variable x
	variable G

	set x []
	set G []

	set pmf [open $pmfFile r]

	while {[gets $pmf line] >= 0} {
		
		if {$line == ""} {
			continue
		}

		set splitedline [splitline $line]

		# the comments of NAMD pmf file start with a "#"
		if {[lindex $line 0] == "#"} {
			continue
		}

		lappend x [lindex $line 0]
		lappend G [lindex $line 1]
	}
}

# find the minimun of pmf
proc ::BFEE::findminimum {} {

	variable x
	variable G
	variable width

	# find the minimun of pmf
	set minx [lindex $x 0]
	set miny [lindex $G 0]
	foreach i $x j $G {
		if {$j < $miny} {
			set miny $j
			set minx $i
		}
	}
	return $minx
}

# do integration, calculate the contribution of RMSD, euler angles and polar angles
proc ::BFEE::integration {k rmsd unbound} {

	variable x
	variable G
	variable width

	# find the minimun of pmf
	if {$rmsd != 1} {
		set minx [lindex $x 0]
		set miny [lindex $G 0]
		foreach i $x j $G {
			if {$j < $miny} {
				set miny $j
				set minx $i
			}
		}
	} else {
		set minx 0
	}

	# integration

	set width [expr [lindex $x 1] - [lindex $x 0]]
	set numerator 0
	set denominator 0
	foreach i $x j $G {
		set numerator [expr $numerator + exp((-[beta]) * $j)]
		set denominator [expr $denominator + exp((-[beta]) * ($j + 0.5 * $k * (($i - $minx)**2)))]
	}

	if {$unbound == 1} {
		return [expr log($numerator / $denominator) / [beta]]
	}
	return [expr -log($numerator / $denominator) / [beta]]
}

# calculate the contribution of euler angles in bulk
proc ::BFEE::calcRigidBody {eulerTheta0 eulerPhi0 eulerPsi0 kTheta kPhi kPsi} {
	
	set eulerTheta0 [expr $eulerTheta0 * 3.1415926 / 180]
	#set eulerPhi0 [expr $eulerPhi0 * 3.1415926 / 180]
	#set eulerPsi0 [expr $eulerPsi0 * 3.1415926 / 180]

	# colvars should have pbc, then the u(phi) and u(psi) should be the same in all cases
	set eulerPhi0 [expr 180 * 3.1415926 / 180]
	set eulerPsi0 [expr 180 * 3.1415926 / 180]
	set kTheta [expr $kTheta * (180/3.1415926)**2]
	set kPhi [expr $kPhi * (180/3.1415926)**2]
	set kPsi [expr $kPsi * (180/3.1415926)**2]
		
	# integration of Theta
	set contriTheta 0
	for {set i -1.570} {$i < 1.570} {set i [expr $i + 0.001]} {
		set contriTheta [expr $contriTheta + 0.001 * sin($i+1.570) * exp(-[beta] * 0.5 * $kTheta * ($i - $eulerTheta0)**2)]
	}

	# integration of Phi
	# these two integration should be irrelevant to the PMF!
	set contriPhi 0
	for {set i 0} {$i < 6.283} {set i [expr $i + 0.001]} {
		set contriPhi [expr $contriPhi + 0.001 * exp(-[beta] * 0.5 * $kPhi * ($i - $eulerPhi0)**2)]
	}

	# integration of Psi
	set contriPsi 0
	for {set i 0} {$i < 6.283} {set i [expr $i + 0.001]} {
		set contriPsi [expr $contriPsi + 0.001 * exp(-[beta] * 0.5 * $kPsi * ($i - $eulerPsi0)**2)]
	}

	return [expr -log(($contriTheta * $contriPhi * $contriPsi) / 8 / (3.1415926**2)) / [beta]]
}
	
# Jacobian correction of the seperation PMF
proc ::BFEE::correctJ {} {
		
	variable x
	variable G
	variable BOLTZMANN
	variable analyzeTemperature
		
	set Gcorr []
	foreach i $x j $G {
		lappend Gcorr [expr $j + 2 * $BOLTZMANN * $analyzeTemperature * log($i)]
	}

	set G $Gcorr
}

# calculate the contribution of S* and I*
proc ::BFEE::calcSI {r polarTheta0 polarPhi0 kTheta kPhi} {
		
	variable x
	variable G

	# if r* > r(max), raise a warning and set r* = rmax
	if {$r > [lindex $x end]} {
		tk_messageBox -type ok -icon info -message \
"Warning: r* > r(max), please check your settings!
Now set r* = r(max) = [lindex $x end] and continue calculation" 
	set r [lindex $x end]
	}

	set polarTheta0 [expr $polarTheta0 * 3.1415926 / 180]
	#set polarPhi0 [expr $polarPhi0 * 3.1415926 / 180]

	# colvars should have pbc, then the u(Phi) should be the same in all cases
	set polarPhi0 [expr 0 * 3.1415926 / 180]
	set kTheta [expr $kTheta * (180/3.1415926)**2]
	set kPhi [expr $kPhi * (180/3.1415926)**2]

	# integration of Theta
	set contriTheta 0
	for {set i 0} {$i < 3.141} {set i [expr $i + 0.001]} {
		set contriTheta [expr $contriTheta + 0.001 * sin($i) * exp(-[beta] * 0.5 * $kTheta * ($i - $polarTheta0)**2)]
	}

	# integration of Phi
	# the colvar difinition of polarPhi is [-3.141, 3.141]
	# due to the periodicity, this value should be irrelevant to the PMF
	set contriPhi 0
	for {set i -3.141} {$i < 3.141} {set i [expr $i + 0.001]} {
		set contriPhi [expr $contriPhi + 0.001 * exp(-[beta] * 0.5 * $kPhi * ($i - $polarPhi0)**2)]
	}

	set S [expr $r * $r * $contriTheta * $contriPhi]

	# w(r*)
	set Gr [lindex $G 0]
	foreach i $x j $G {
		if {$i >= $r} {
			set Gr $j
		}
	}
	# integration of r
	set width [expr [lindex $x 1] - [lindex $x 0]]
	set I 0
	foreach i $x j $G {
		set I [expr $I + $width * exp(-[beta] * ($j - $Gr))]
	}

	return [expr -1 / [beta] * log($S * $I / 1661)]
}

# calculate binding free energy
proc ::BFEE::calculate_free_energy {dir1 fc1 dir2 fc2 dir3 fc3 dir4 fc4 dir5 fc5 dir6 fc6
dir7 r_star dir8 fc8} {

	readPMF $dir1
	set dg1 [integration $fc1 1 0]

	readPMF $dir2
	set dg2 [integration $fc2 0 0]
	set minx2 [findminimum]

	readPMF $dir3
	set dg3 [integration $fc3 0 0]
	set minx3 [findminimum]

	readPMF $dir4
	set dg4 [integration $fc4 0 0]
	set minx4 [findminimum]

	readPMF $dir5
	set dg5 [integration $fc5 0 0]
	set minx5 [findminimum]

	readPMF $dir6
	set dg6 [integration $fc6 0 0]
	set minx6 [findminimum]

	readPMF $dir8
	set dg8 [integration $fc8 1 1]

	readPMF $dir7
	correctJ
	set dg7 [calcSI $r_star $minx5 $minx6 $fc5 $fc6]

	tk_messageBox -type ok -icon info \
	-message \
"[format "%-25s" dG(site,c)] = [format "%+10.2f" $dg1] kcal/mol
[format "%-25s" dG(site,eulerTheta)] = [format "%+10.2f" $dg2] kcal/mol
[format "%-25s" dG(site,eulerPhi)] = [format "%+10.2f" $dg3] kcal/mol
[format "%-25s" dG(site,eulerPsi)] = [format "%+10.2f" $dg4] kcal/mol
[format "%-25s" dG(site,polarTheta)] = [format "%+10.2f" $dg5] kcal/mol
[format "%-25s" dG(site,polarPhi)] = [format "%+10.2f" $dg6] kcal/mol
[format "%-25s" (1/beta)*ln(S*I*C0)] = [format "%+10.2f" $dg7] kcal/mol
[format "%-25s" dG(bulk,c)] = [format "%+10.2f" $dg8] kcal/mol
[format "%-25s" dG(bulk,o)] = [format "%+10.2f" [set dg9 [calcRigidBody $minx2 $minx3 $minx4 $fc2 $fc3 $fc4]]] kcal/mol

Standard Binding Free Energy:
[format "%-25s" dG(total)] = [format "%+.2f" [expr $dg1 + $dg2 + $dg3 + $dg4 + $dg5 + $dg6 + $dg7 + $dg8 + $dg9]] kcal/mol"
}

# compute Binding Free energy
proc ::BFEE::computeFreeEnergy {} {
	variable rmsdDir
	variable analyzermse 
	variable eulerThetaDir 
	variable analyzeeulerTheta 
	variable eulerPhiDir 
	variable analyzeeulerPhi 
	variable eulerPsiDir 
	variable analyzeeulerPsi
	variable polarThetaDir 
	variable analyzepolarTheta 
	variable polarPhiDir 
	variable analyzepolarPhi 
	variable rDir
	variable unrmsdDir
	variable analyzermsd
	variable r_star

	return [calculate_free_energy $rmsdDir $analyzermsd $eulerThetaDir $analyzeeulerTheta $eulerPhiDir $analyzeeulerPhi $eulerPsiDir $analyzeeulerPsi \
	$polarThetaDir $analyzepolarTheta $polarPhiDir $analyzepolarPhi $rDir $r_star $unrmsdDir $analyzermsd]
}
}

proc bfee_tk {} {
	::BFEE::drawGUI
	return $::BFEE::w
}