# Copyright (c) 2010, University of Pittsburgh
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIlow. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

package provide comd 1.0

package require solvate
package require autoionize
package require psfgen
package require pbctools
package require exectool


set COMD_PATH $env(COMD_PATH)
set PACKAGE_PATH "$COMD_PATH"
set PACKAGEPATH "$COMD_PATH"
set env(PROBEDIR) "$COMD_PATH"

variable platform $tcl_platform(platform)
switch $platform {
  unix {
    set TMPDIR "/tmp" ;  # or even $::env(TMPDIR), at times.
  } macintosh {
    set TMPDIR $::env(TRASH_FOLDER)  ;# a better place?
  } default {
    set TMPDIR [pwd]
    catch {set TMPDIR $::env(TMP)}
    catch {set TMPDIR $::env(TEMP)}
  }
}

namespace eval ::comd:: {
  namespace export comd

  variable version 1.0

  variable w

  # Variables for system setup
  variable molid -1
  # input files
  variable initial_pdb 
  variable final_pdb 
  variable initial_chid 
  variable final_chid 
  # Ionization parameters
  variable topo_file 
  variable solvent_padding 
  variable neutralize 
  # Minimization parameters
  variable para_file 
  variable temperature 
  variable min_length 
  # ANM-MC-Metropolis parameters
  variable anm_cutoff 
  variable dev_mag 
  # TMD options
  variable spring_k 
  variable tmd_len 
  # Simulation options
  variable comd_cycle 
  variable num_cores 
  variable python_path ""
  # output options
  variable outputdir 
  variable output_prefix 
  
  # Logvew window counter
  variable logcount 0
  variable lognames [list]
  variable titles [list "Prepare System"]
  variable interfaces [list "prepare"]
  variable which_mode [lindex $titles 0]


}


proc comd::Logview {logfilename} {
  variable logcount
  variable lognames
  set logindex [lsearch $lognames $logfilename]
  set log .somenonsense
  if {$logindex > -1} {
    set windowname "log$logindex"
    set log .$windowname
  }
  if {[winfo exists $log] == 0} {
    if {$logindex > -1} {
      lset lognames $logindex "somenonsense"
    }
    set logindex $logcount
    lappend lognames $logfilename
    set windowname "log$logindex"
    set log [toplevel ".$windowname"]
    wm title $log "Logfile [lindex [file split $logfilename] end] ($logfilename)"
    wm resizable $log 1 1
    incr logcount

    text $log.text -bg White -bd 2 \
      -yscrollcommand ".$windowname.vscr set"
    scrollbar $log.vscr -command ".$windowname.text yview"
    pack $log.text -side left -fill both -expand 1
    pack $log.vscr -side right -fill y
  }

  $log.text configure -state normal
  #set count 0
  #set tabwidth 0
  #foreach family [lsort -dictionary [font families]] {
  #    $log.text tag configure f[incr count] -font [list $family 10]
  #    $log.text insert end ${family}:\t {} \
  #            "This is a simple sampler\n" f$count
  #    set w [font measure [$log.text cget -font] ${family}:]
  #    if {$w+5 > $tabwidth} {
  #        set tabwidth [expr {$w+5}]
  #        $log.text configure -tabs $tabwidth
  #    }
  #}
  $log.text delete 1.0 end
  set logfile [open $logfilename "r"]
  set line ""
  while {[gets $logfile line] != -1} {
    $log.text insert end "$line\n"
  }
  close $logfile
  $log.text yview moveto 1
  $log.text configure -state disabled
}

proc ::comd::comdgui {} {
  variable w

  global env

  # Determine whether PSF and PDB are loaded (based on Solvate plugin code)
  # set ::comd::initial_pdb ""
  # set ::comd::final_pdb ""


  # If already initialized, just turn on
  if [winfo exists .comdgui] {
    wm deiconify .comdgui
    raise .comdgui
    return
  }

  # Initialize window
  set w [toplevel .comdgui]
  wm title $w "COllective Molecular Dynamics v$::comd::version"
  wm resizable $w 0 0

#   set wif [frame $w.interface_frame]
#   button $wif.help -text "?" -padx 0 -pady 3 -command {
#       tk_messageBox -type ok -title "HELP" \
#       -message "Use option menu to change the active interface. There are\
# two interfaces to facilitate the COllective Molecular Dynamics.\n\n\
# [lindex $::druggability::titles 0]\n\
# Prepare protein in a water-probe mixture box or in a water-only box using\
# this interface. Also, by default generic NAMD input files are outputed.\
# Protein PSF and PDB files that also \
# contains cofactors/ions etc. are required from the user.\n\n\
# "}
#   variable titles
#   tk_optionMenu $wif.list ::comd::which_mode "System Setup"
#   $wif.list.menu delete 0
#   $wif.list.menu add radiobutton -label [lindex $titles 0] \
#     -variable ::comd::which_mode \
#     -command {::comd::Switch_mode "prepare"}
#   pack $wif.help -side left
#   pack $wif.list -side left -expand 1 -fill x
#   pack $wif -pady 2 -expand 1 -fill x

  # Set main frame
  set mf [frame $w.main_frame]

  # VISUALIZE results
  
  # Prepare System and Simulation Files
  set mfa [frame $mf.prepare]
  # Select input files
  set mfaif [labelframe $mfa.input_files -text "Protein structures:" -bd 2]
  # Initial PDB
  grid [button $mfaif.psf_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The initial protein structure should be given in the standard PDB format."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaif.psf_label -text "Initial PDB: "] -row 1 -column 1 -sticky w
  grid [entry $mfaif.psf_path -width 28 \
      -textvariable ::comd::initial_pdb] \
    -row 1 -column 2 -columnspan 4 -sticky ew
  grid [button $mfaif.psf_browse -text "Browse" -width 6 -pady 1 -command {
      set tempfile [tk_getOpenFile \
                    -filetypes {{"PDB files" { .pdb .PDB }} {"All files" *}}]
      if {![string equal $tempfile ""]} {
        set ::comd::initial_pdb $tempfile
      } }] \
    -row 1 -column 6 -sticky w
       
  # Final PDB
  grid [button $mfaif.pdb_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The final protein structure should be given in the standard PDB format."}] \
    -row 2 -column 0 -sticky w
  grid [label $mfaif.pdb_label -text "Final PDB: "] \
    -row 2 -column 1 -sticky w
  grid [entry $mfaif.pdb_path -width 28 \
      -textvariable ::comd::final_pdb] \
    -row 2 -column 2 -columnspan 4 -sticky ew
  grid [button $mfaif.pdb_browse -text "Browse" -width 6 -pady 1 -command {
        set tempfile [tk_getOpenFile \
          -filetypes {{"PDB files" { .pdb .PDB }} {"All files" *}}]
        if {![string equal $tempfile ""]} {
          set ::comd::final_pdb $tempfile
        } }] \
    -row 2 -column 6 -sticky w

  grid [button $mfaif.inich_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The chain ID for the initial and final structure which should be in the previously imported PDB file."}] \
    -row 3 -column 0 -sticky w
  grid [label $mfaif.inich_label -text "Initial PDB chain ID: "] \
    -row 3 -column 1 -sticky w
  grid [entry $mfaif.inich_entry -width 6 \
    -textvariable ::comd::initial_chid] \
    -row 3 -column 2 -sticky ew  

  #grid [button $mfaif.finch_help -text "?" -padx 0 -pady 0 -command {
  #    tk_messageBox -type ok -title "HELP" \
        -message "The chain ID for the initial structure which should be in the previously imported PDB file."}] \
    -row 3 -column 4 -sticky w
  grid [label $mfaif.finch_label -text "Final PDB chain ID: "] \
    -row 3 -column 4 -columnspan 2 -sticky w
  grid [entry $mfaif.finch_entry -width 6 \
    -textvariable ::comd::final_chid] \
    -row 3 -column 6 -sticky ew
  
    
  pack $mfaif -side top -ipadx 0 -ipady 5 -fill x -expand 1

  # Enter ionization options
  set mfaio [labelframe $mfa.ionize_options -text "Ionization parameters" -bd 2]
  # Topology File
  grid [button $mfaio.topo_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "This file is the topology file that will be given to psfgen package of vmd. Therefore, \
        from given static structure this will create bonds, angles and various structural elements based on \
        topology parameters provided. Suggested file extension is .top but others will be accepted. "}] \
    -row 0 -column 0 -sticky w
  grid [label $mfaio.topo_label -text "Topology File:"] -row 0 -column 1 -sticky w
  grid [entry $mfaio.topo_path -width 28 \
      -textvariable ::comd::topo_file] \
    -row 0 -column 2 -columnspan 4 -sticky ew
  grid [button $mfaio.topo_browse -text "Browse" -width 6 -pady 1 -command {
      set tempfile [tk_getOpenFile \
                    -filetypes {{"Topology files" { .top .TOP }} {"All files" *}}]
      if {![string equal $tempfile ""]} {
        set ::comd::topo_file $tempfile
      } }] \
    -row 0 -column 6 -sticky w
  
  #Solvation box padding and counter ions
  grid [button $mfaio.padding_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "This is the half of the initial distance between the protein\
and its imaginary copies under periodic boundary conditions. For systems with\
probes, the resulting padding distance will be slightly larger, due to\
constraint of preserving the ratio of 20 water molecules per probe molecule."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaio.padding_label -text "Box padding (A): "] \
    -row 1 -column 1 -sticky w
  grid [entry $mfaio.padding_entry -width 6 \
    -textvariable ::comd::solvent_padding] \
    -row 1 -column 2 -sticky ew

  #grid [label $mfaio.separatpr_label -text "   "] \
    -row 1 -column 3 -sticky w

  grid [button $mfaio.neutralize_help -text "?" -padx 0 -pady 0 -command {
    tk_messageBox -type ok -title "HELP" \
      -message "By default, counter ions will be added to neutralize a charged\
system. A charged system (if the protein is charged) may be obtained by unchecking this option."}] \
    -row 1 -column 4 -sticky w
  grid [label $mfaio.neutralize_label \
      -text "Add counter ions: "] \
    -row 1 -column 5 -sticky w
  grid [checkbutton $mfaio.neutralize_check -text "" \
      -variable ::comd::neutralize] \
    -row 1 -column 6 -sticky w

  pack $mfaio -side top -ipadx 0 -ipady 5 -fill x -expand 1

  # Enter minimization options
  set mfamo [labelframe $mfa.minimize_options -text "Minimization parameters" -bd 2]
  # Topology File
  grid [button $mfamo.para_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The parameter files for force field. The file should be provided in par or prm format and include necessary parameters required for NAMD."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfamo.para_label -text "Parameter File: "] -row 0 -column 1 -sticky w
  grid [entry $mfamo.para_path -width 28 \
      -textvariable ::comd::para_file] \
    -row 0 -column 2 -columnspan 4 -sticky ew
  grid [button $mfamo.para_browse -text "Browse" -width 6 -pady 1 -command {
      set tempfile [tk_getOpenFile \
                    -filetypes {{"Parameter files" { .par .PAR .prm .PRM }} {"All files" *}}]
      if {![string equal $tempfile ""]} {
        set ::comd::para_file $tempfile
      } }] \
    -row 0 -column 6 -sticky w
  
  #Temperature and minimization length parameters
  grid [button $mfamo.temperature_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The temperature for molecular dynamics simulation needs to be entered. The units are in Kelvin."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfamo.temperature_label -text "Temperature (K): "] \
    -row 1 -column 1 -sticky w
  grid [entry $mfamo.temperature_entry -width 6 \
    -textvariable ::comd::temperature] \
    -row 1 -column 2 -sticky ew

  grid [label $mfamo.separatpr_label -text "      "] \
    -row 1 -column 3 -sticky w

  grid [button $mfamo.length_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "After running targeted molecular dynamics simulations the final structure needs to be equilibrated into a stable state. The length of minimization creates 
        more stable structures which will guarantees user to have a stable targeted molecular dynamics simulation. The units are in ps. "}] \
    -row 1 -column 4 -sticky w
  grid [label $mfamo.length_label -text "Minimization length (ps): "] \
    -row 1 -column 5 -sticky w
  grid [entry $mfamo.length_entry -width 3 \
    -textvariable ::comd::min_length] \
    -row 1 -column 6 -sticky ew
  

  pack $mfamo -side top -ipadx 0 -ipady 5 -fill x -expand 1

  set mfamc [labelframe $mfa.anmmc_options -text "ANM-MC-Metropolis options:" -bd 2]
  
  grid [button $mfamc.anmcut_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "In ANM calculations, the cutoff parameter is the maximum distance that two proteins are in contact. The units are A. "}] \
    -row 0 -column 0 -sticky w
  grid [label $mfamc.anmc_label -text "ANM cutoff (A): "] \
    -row 0 -column 1 -sticky w
  grid [entry $mfamc.anmc_field -width 6 -textvariable ::comd::anm_cutoff] \
    -row 0 -column 2 -sticky w

  grid [label $mfamc.separatpr2_label -text "      "] \
    -row 0 -column 3 -sticky w

  grid [button $mfamc.dev_mag_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The scaling factor used when disturbing the structure of protein in ANM-MC step. Default and suggested value is 0.1 A."}] \
    -row 0 -column 4 -sticky w
  grid [label $mfamc.dev_mag_label -text "Deviation (A): "] \
    -row 0 -column 5 -sticky w
  grid [entry $mfamc.dev_mag_field -width 6 \
      -textvariable ::comd::dev_mag] \
    -row 0 -column 6 -sticky w



  # grid [button $mfamc.stepcut_help -text "?" -padx 0 -pady 0 -command {
  #     tk_messageBox -type ok -title "HELP" \
  #       -message "To keep structure intact and to avoid having unrealistic and very different structures in ANM-MC step, an rmsd threshold is used. Suggested value is 4 A."}] \
  #   -row 1 -column 0 -sticky w
  # grid [label $mfamc.stepc_label -text "Step cutoff (A): "] \
  #   -row 1 -column 1 -sticky w
  # grid [entry $mfamc.stepc_field -width 6 -textvariable ::comd::step_cutoff] \
  #   -row 1 -column 2 -sticky w

  #   grid [label $mfamc.separatpr_label -text "      "] \
  #   -row 1 -column 3 -sticky w

  # grid [button $mfamc.num_cyc_help -text "?" -padx 0 -pady 0 -command {
  #     tk_messageBox -type ok -title "HELP" \
  #       -message "The number of disturbance on the structure based on ANM modes. The number of disturbance should be selected based on structure and for a structure with 200 residues suggested number is in order of thousands."}] \
  #   -row 1 -column 4 -sticky w
  # grid [label $mfamc.num_cyc_label -text "No of disturbance: "] \
  #   -row 1 -column 5 -sticky w
  # grid [entry $mfamc.num_cyc_field -width 6 \
  #     -textvariable ::comd::anm_cycle] \
  #   -row 1 -column 6 -sticky w  

  pack $mfamc -side top -ipadx 0 -ipady 5 -fill x -expand 1

  set mfatm [labelframe $mfa.tmd_options -text "TMD options:" -bd 2]

  grid [button $mfatm.spring_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "In targeted molecular dynamics simulation, the target potential is harmonic and the spring constant term shows the force applied to structure to reach final structure."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfatm.spring_label -text "Spring constant:"] \
    -row 0 -column 1 -sticky w
  grid [entry $mfatm.spring_field -width 6 -textvariable ::comd::spring_k] \
    -row 0 -column 2 -sticky w

  grid [label $mfatm.separatpr2_label -text "   "] \
    -row 0 -column 3 -sticky w

  grid [button $mfatm.tmd_len_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The length of targeted molecular dynamics simulations in the units of ps. The length of collective molecular dynamics will change based on the structure and for a structure with 200 residues suggested length is in the order of hundreds."}] \
    -row 0 -column 4 -sticky w
  grid [label $mfatm.tmd_len_label -text "TMD length:"] \
    -row 0 -column 5 -sticky w
  grid [entry $mfatm.tmd_len_field -width 6 \
      -textvariable ::comd::tmd_len] \
    -row 0 -column 6 -sticky w
  
  pack $mfatm -side top -ipadx 0 -ipady 5 -fill x -expand 1

  ########################################3
  set mfaso [labelframe $mfa.simulation_options -text "Simulation options:" -bd 2]

  grid [button $mfaso.cmd_cyc_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "Each CoMD cycle consists of minimization, ANM-MC-Metropolis distrubance and targeted molecular dynamics. The total number of cycle gives the maximum number of cycles perforlow before starting and final structures are very close."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfaso.cmd_cyc_label -text "No of coMD cycles:"] \
    -row 0 -column 1 -sticky w
  grid [entry $mfaso.cmd_cyc_field -width 6 -textvariable ::comd::comd_cycle] \
    -row 0 -column 2 -sticky w

  grid [label $mfaso.separatpr2_label -text "   "] \
    -row 0 -column 3 -sticky w

  grid [button $mfaso.numcores_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "The number of physical cores in the cluster or PC that you will run your simulation. NAMD is running parallel on CPUs."}] \
    -row 0 -column 3 -sticky w
  grid [label $mfaso.num_cores_label -text "No of physical cores:"] \
    -row 0 -column 4 -sticky w
  grid [entry $mfaso.num_cores -width 5 \
      -textvariable ::comd::num_cores] \
    -row 0 -column 5 -columnspan 4 -sticky ew
  
  pack $mfaso -side top -ipadx 0 -ipady 5 -fill x -expand 1



  ########################################4


  set mfaoo [labelframe $mfa.output_options -text "Output options:" -bd 2]

  grid [button $mfaoo.outdir_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "Output folder, default is current working directory."}] \
    -row 0 -column 0 -sticky w
  grid [label $mfaoo.outdir_label -text "Output folder:"] \
    -row 0 -column 1 -sticky w
  grid [entry $mfaoo.outdir_path -width 28 \
      -textvariable ::comd::outputdir] \
    -row 0 -column 2 -columnspan 4 -sticky ew
  grid [button $mfaoo.dcd_browse -text "Browse" -width 6 -pady 1 -command {
      set tempfile [tk_chooseDirectory]
      if {![string equal $tempfile ""]} {
        set ::comd::outputdir $tempfile
      }}] \
    -row 0 -column 6 -sticky w

  grid [button $mfaoo.prefix_help -text "?" -padx 0 -pady 0 -command {
      tk_messageBox -type ok -title "HELP" \
        -message "All output files and folders will start with this prefix.\
A unique and descriptive prefix choice may allow running multiple simulations in the same folder."}] \
    -row 1 -column 0 -sticky w
  grid [label $mfaoo.prefix_label -text "Output prefix:"] \
    -row 1 -column 1 -sticky w
  grid [entry $mfaoo.prefix_path -width 28 \
      -textvariable ::comd::output_prefix] \
    -row 1 -column 2 -sticky w

  grid [label $mfaoo.separator_label -text "   "] \
    -row 1 -column 3 -columnspan 4 -sticky w
  
  pack $mfaoo -side top -ipadx 0 -ipady 5 -fill x -expand 1

  # Prepare System
  button $mfa.button -text "Prepare System" -command ::comd::Prepare_system -bd 3
  pack $mfa.button

  pack $mfa -side top -padx 0 -pady 0 -fill x -expand 1
  pack $mf -side top -padx 0 -pady 0 -fill x -expand 1

  return $w
}


# proc ::druggability::Switch_mode {which_mode} {
#   # change GUI layout
#   variable w
#   pack forget $w.main_frame.prepare $w.main_frame.process $w.main_frame.analyze $w.main_frame.evaluate $w.main_frame.visualize
#   pack $w.main_frame.$which_mode -side top -padx 0 -pady 0 -fill x -expand 1
#   variable interface $which_mode
# }

# proc ::druggability::Load_protein {} {
#   # load protein and set molid
#   variable initial_pro
#   variable final_pro
#   variable initial_chid
#   variable final_chid
#   if {$::druggability::initial_pdb != "" && $::druggability::final_pdb != ""} {
#     mol new $::druggability::initial_pdb 
#     set initial_pro [atomselect top "protein and not altloc B and not hydrogen and chain $initial_chid"]
#     mol new $::druggability::final_pdb 
#     set final_pro [atomselect top "protein and not altloc B and not hydrogen and chain $final_chid"]
#     return 1
#   } else {
#     tk_messageBox -type ok -title "ERROR" \
#       -message "Both PDB files must be specified."
#     return 0
#   }
# }

proc ::comd::Prepare_system {} {

  # WHAT IS NEW?
  # 2.1 - Bug fixes, and file checks
  # 2.0 - Allows setup of systems containing multiple probe tybes
  # 2.0 - Improved system setup provides lesser number of solvent atoms
  # 2.0 - Cleans up interlowiate files
  # 2.0 - Outputs a log file for trouble shooting, and further intstructions
  # 2.0 - NAMD configuration files are prepared for a single or multiple
  #       simulations
  # HOW THE CODE WORKS
  # The code will
  # (1)   solvate the protein, or everything in the PDB/PSF files that you provide
  # (2)   add counter ions to neutralize the system
  # (3)   adjust the water/probe ratio to a predefined level (1 probe to 20 water)
  # (4)   write output files for each simula tion
  #       Setups of multiple simulations differ only at random number seeds.
  #       This will be sufficient to result in a different trajectory.
  variable w
  variable initial_pdb
  variable final_pdb
  variable initial_chid
  variable final_chid
  variable pro
  variable topo_file
  variable temperature
  variable para_file
  variable min_length
  variable percent_ipro
  variable percent_ibut
  variable percent_acam
  variable percent_acetipam
  variable solvent_padding
  variable neutralize
  variable output_prefix
  variable comd_cycle
  variable anm_cutoff 
  variable dev_mag 
  variable spring_k
  variable tmd_len
  variable num_cores
  variable outputdir
  variable python_path

  if {$outputdir != ""} {
      if {![file isdirectory $outputdir]} {
        if {[catch {file mkdir $outputdir}]} {
          tk_messageBox -type ok -title "ERROR" \
            -message "Could not make output folder: $outputdir"
          return

        }
      }
  }

  if {$::comd::initial_pdb == "" || $::comd::final_pdb == ""} {
    tk_messageBox -type ok -title "ERROR" \
      -message "Both PDB files must be specified."
    return
  }

  #set percent_total [expr $percent_acam + $percent_ibut + $percent_ipro + $percent_acetipam]

  if {$solvent_padding < 4} {
    tk_messageBox -type ok -title "ERROR" \
      -message "Solvent box padding parameter must be larger than 4 A."
    return
  }

  if {[string length [string trim $output_prefix]] == 0} {
    tk_messageBox -type ok -title "ERROR" \
      -message "Please enter a descriptive name (prefix) for the system."
    return
  }

  global env
  global COMD_PATH

  set log_file [open [file join "$outputdir" "$output_prefix.log"] a]
  puts $log_file "---==## [clock format [clock seconds]] #==---"
  puts $log_file "Version: $::comd::version"
  puts $log_file "Info: Logging started for setup of $output_prefix."
  puts $log_file "Solvation: Box padding $solvent_padding A."
  
  package require psfgen
  package require solvate
  package require autoionize

  
  ####### SOLVATION AND IONIZATION OF INITIAL PROTEIN STRUCTURE #######


  resetpsf
  mol delete all
  mol new $initial_pdb
  set pro [atomselect top "protein and chain ${initial_chid}"]
  $pro writepdb init.pdb
  mol delete all
  mol new init.pdb
  set pro [atomselect top "protein and not altloc B and not hydrogen and chain ${initial_chid}"]
  $pro writepdb prop.pdb

  if {$topo_file == ""} {
  set topo_file "${COMD_PATH}/top_all27_prot_lipid.top"
  }
  topology $topo_file
  pdbalias residue HIS HSD
  pdbalias atom ILE CD1 CD
  segment R {pdb prop.pdb}
  coordpdb prop.pdb R
  guesscoord
  writepdb pro.pdb
  writepsf pro.psf
  solvate pro.psf pro.pdb -t $solvent_padding -o pro_wb
  # DELETE solvated molecule

  if {$neutralize} {
    set totalcharge 0
    foreach charge [[atomselect top "all"] get charge] {
      set totalcharge [expr $totalcharge + $charge]
    }
    # number of CL and NA atoms are determined
    puts $log_file "Ionization: Initial PDB System has a total charge of $totalcharge electrons."
    if {$totalcharge > 0} {
        set nna 0
        set ncl [expr round($totalcharge)]
        puts $log_file "Ionization: $ncl chloride ions will be added to initial PDB."
    } else {
        set ncl 0
        set nna [expr -1 * round($totalcharge)]
        puts $log_file "Ionization: $nna sodium ions will be added to initial PDB."
    }
    if {$ncl > 0 | $nna > 0} {
        autoionize -psf pro_wb.psf -pdb pro_wb.pdb \
        -o [file join ${outputdir} "initial_ionized"] -from 5 -between 5 -ncl $ncl -nna $nna -seg ION
        puts $log_file "Ionization: Initial PDB System is ionized to become neutral."
    }
  }
  
  mol new init.pdb
  set everyone [atomselect top "all"]
  set initial_lims [measure minmax $everyone]
  set initial_cent [measure center $everyone]
  set ixmin [lindex $initial_lims 0 0]
  set ixmax [lindex $initial_lims 1 0]
  set iymin [lindex $initial_lims 0 1]
  set iymax [lindex $initial_lims 1 1]
  set izmin [lindex $initial_lims 0 2]
  set izmax [lindex $initial_lims 1 2]
  set ixcen [lindex $initial_cent 0]
  set iycen [lindex $initial_cent 1]
  set izcen [lindex $initial_cent 2]
  set ixlen [expr {$ixmax-$ixmin+16}]
  set iylen [expr {$iymax-$iymin+16}]
  set izlen [expr {$izmax-$izmin+16}]


    ####### SOLVATION AND IONIZATION OF FINAL PROTEIN STRUCTURE #######

  resetpsf
  mol delete all
  mol new $final_pdb
  set pro [atomselect top "protein and chain ${final_chid}"]
  $pro writepdb fino.pdb
  mol delete all
  mol new fino.pdb
  set pro [atomselect top "protein and not altloc B and not hydrogen and chain ${final_chid}"]
  $pro writepdb prop.pdb
  topology $topo_file
  pdbalias residue HIS HSD
  pdbalias atom ILE CD1 CD1
  segment R {pdb prop.pdb}
  coordpdb prop.pdb R
  guesscoord
  writepdb pro.pdb
  writepsf pro.psf
  solvate pro.psf pro.pdb -t $solvent_padding -o pro_wb

  # # DELETE solvated molecule

  if {$neutralize} {
    set totalcharge 0
    foreach charge [[atomselect top "all"] get charge] {
      set totalcharge [expr $totalcharge + $charge]
    }
    # number of CL and NA atoms are determined
    puts $log_file "Ionization: Final PDB System has a total charge of $totalcharge electrons."
    if {$totalcharge > 0} {
        set nna 0
        set ncl [expr round($totalcharge)]
        puts $log_file "Ionization: $ncl chloride ions will be added to final PDB."
    } else {
        set ncl 0
        set nna [expr -1 * round($totalcharge)]
        puts $log_file "Ionization: $nna sodium ions will be added to final PDB."
    }
    if {$ncl > 0 | $nna > 0} {
        autoionize -psf pro_wb.psf -pdb pro_wb.pdb \
        -o [file join ${outputdir} "final_ionized"] -from 5 -between 5 -ncl $ncl -nna $nna -seg ION
        puts $log_file "Ionization: Final PDB System is ionized to become neutral."
    }
  }

  mol new fino.pdb
  set everyone [atomselect top "all"]
  set final_lims [measure minmax $everyone]
  set final_cent [measure center $everyone]
  set fxmin [lindex $final_lims 0 0]
  set fxmax [lindex $final_lims 1 0]
  set fymin [lindex $final_lims 0 1]
  set fymax [lindex $final_lims 1 1]
  set fzmin [lindex $final_lims 0 2]
  set fzmax [lindex $final_lims 1 2]
  set fxcen [lindex $final_cent 0]
  set fycen [lindex $final_cent 1]
  set fzcen [lindex $final_cent 2]
  set fxlen [expr {$fxmax-$fxmin+16.0}]
  set fylen [expr {$fymax-$fymin+16.0}]
  set fzlen [expr {$fzmax-$fzmin+16.0}]
  
    ####### INITIAL MINIMIZATION OF INITIAL PROTEIN STRUCTURE #######

  # Initial minimization
  #set status [exec bash "${sh_filename}"]
  #puts $log_file "Simulation: NAMD configuration files for minimization are written into folder $output_prefix$minfix."

  if {$para_file == ""} {
    set para_file "${COMD_PATH}/par_all27_prot_lipid.prm"
  }
  set tcl_file [open [file join "$outputdir" "$output_prefix.tcl"] w]
  #set tcl_file_name [file join "$outputdir" "$output_prefix.tcl"]
  puts $tcl_file "#This tcl file will run full collective molecular dynamics simulation with given parameters."

  #### #
  puts $tcl_file "set sh_file \[open \"$output_prefix.sh\" w\]"
  puts $tcl_file "set sh_filename \"${output_prefix}.sh\""
  puts $tcl_file "package require exectool"
  puts $tcl_file "set namd2path \[::ExecTool::find \"namd2\"\]"
  if {$::comd::python_path == ""} {
  	puts $tcl_file "set python_path \[::ExecTool::find \"python\"\]"	
  } else {
  	puts $tcl_file "set python_path ${python_path}\/python" 
  }
  puts $tcl_file "puts \$sh_file \"\\\#\\\!\\\/bin\\\/bash\""
  puts $tcl_file "puts \$sh_file \"NAMD=\\\"\$namd2path \+p[expr ${num_cores}/2]\\\"\""
  puts $tcl_file "file mkdir \"${output_prefix}_inimin\""
  puts $tcl_file "set namd_file \[open \[file join \"${output_prefix}_inimin\" \"min.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ..\/initial_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ..\/initial_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname initial_minimized0\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  puts $tcl_file "puts \$namd_file \"parameters ../parameter_file.prm\""
  puts $tcl_file "puts \$namd_file \"temperature \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector1 ${ixlen},0,0\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector2 0,${iylen},0\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector3 0,0,${izlen}\""
  puts $tcl_file "puts \$namd_file \"cellOrigin ${ixcen},${iycen},${izcen}\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"outputPressure [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"restartfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"minimize [expr $min_length*500]\""
  puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
  puts $tcl_file "close \$namd_file"
  
  puts $tcl_file "file mkdir \"${output_prefix}_finmin\""
  puts $tcl_file "set namd_file \[open \[file join \"${output_prefix}_finmin\" \"min.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ..\/final_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ..\/final_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname final_minimized0\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  puts $tcl_file "puts \$namd_file \"parameters ../parameter_file.prm\""
  puts $tcl_file "puts \$namd_file \"temperature \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector1 ${fxlen},0,0\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector2 0,${fylen},0\""
  puts $tcl_file "puts \$namd_file \"cellBasisVector3 0,0,${fzlen}\""
  puts $tcl_file "puts \$namd_file \"cellOrigin ${fxcen},${fycen},${fzcen}\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"outputPressure [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"restartfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"minimize [expr $min_length*500]\""
  puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${output_prefix}_inimin\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD min.conf > min.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\"" 
  puts $tcl_file "puts \$sh_file \"cd ${output_prefix}_finmin\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD min.conf > min.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\"" 
  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"
  puts $tcl_file "set status \[catch \{exec mv \/${output_prefix}_inimin\/initial_minimized0.dcd initial_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv initial_trajectory.dcd initr.dcd\} output\]" 
  puts $tcl_file "set status \[catch \{exec mv \/${output_prefix}_finmin\/final_minimized0.dcd final_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv final_trajectory.dcd fintr.dcd\} output\]" 
  
  puts $tcl_file "package require psfgen"
  puts $tcl_file "mol delete all" 
  puts $tcl_file "mol load psf initial_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_inimin/initial_minimized0.coor" 
  puts $tcl_file "set sel1 \[atomselect top \"name CA\"\]" 
  puts $tcl_file "set sel1a \[atomselect top all\]"
  puts $tcl_file "mol load psf final_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_finmin/final_minimized0.coor"  
  puts $tcl_file "set sel2 \[atomselect top \"name CA\"\]" 
  puts $tcl_file "set sel2a \[atomselect top all\]"
  puts $tcl_file "set trans_mat \[measure fit \$sel2 \$sel1\]"
  puts $tcl_file "\$sel2a move \$trans_mat"
  puts $tcl_file "set rmsd \[measure rmsd \$sel2 \$sel1\]"
  puts $tcl_file "set all_rmsd(0) rmsd"
  puts $tcl_file "puts \$rmsd"

  puts $tcl_file "file mkdir ${output_prefix}_inipro"
  puts $tcl_file "file mkdir ${output_prefix}_finpro"
  #loop start
  puts $tcl_file "for {set cycle 0} {\$cycle < ${comd_cycle}} {incr cycle} {"
  puts $tcl_file "mol delete all"
  puts $tcl_file "resetpsf"
  puts $tcl_file "mol load psf initial_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_inimin/initial_minimized\$cycle.coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set s2 \[atomselect top \"all\"\]"
  puts $tcl_file "mol load psf final_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_finmin/final_minimized\$cycle.coor"
  puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
  puts $tcl_file "\$s2 move \$trans_mat"
  puts $tcl_file "\$s1 writepdb starting_initial.pdb"
  puts $tcl_file "mol delete all"
  puts $tcl_file "resetpsf"
  puts $tcl_file "mol load psf final_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_finmin/final_minimized\$cycle.coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set s2 \[atomselect top \"all\"\]"
  puts $tcl_file "mol load psf initial_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_inimin/initial_minimized\$cycle.coor"
  puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
  puts $tcl_file "\$s2 move \$trans_mat"
  puts $tcl_file "\$s1 writepdb starting_final.pdb"
  
  puts $tcl_file "mol delete all"
  puts $tcl_file "resetpsf"
  puts $tcl_file "mol load psf final_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_finmin/final_minimized\$cycle.coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "\$s1 writepdb initial_target.pdb"
  puts $tcl_file "mol delete all"
  puts $tcl_file "resetpsf"
  puts $tcl_file "mol load psf initial_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_inimin/initial_minimized\$cycle.coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "\$s1 writepdb final_target.pdb"
  
  #set anmmc_path [file join "$COMD_PATH" "anmmc.py"]
  puts $tcl_file "set result \[exec -ignorestderr \$python_path anmmc.py starting_initial.pdb initial_target.pdb ${anm_cutoff} ${dev_mag}\]"
  puts $tcl_file "set result \[exec -ignorestderr \$python_path anmmc.py starting_final.pdb final_target.pdb ${anm_cutoff} ${dev_mag}\]"
  
  puts $tcl_file "mol delete all"
  puts $tcl_file "mol load psf initial_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_inimin/initial_minimized\$cycle.coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set s2 \[atomselect top \"all\"\]"
  puts $tcl_file "mol load pdb final_target.pdb"
  puts $tcl_file "mol addfile starting_initial.pdb_initial_target.pdb_final_structure.dcd"
  puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
  puts $tcl_file "\$s3 move \$trans_mat"
  puts $tcl_file "\$s1 set \{x y z\} \[\$s3 get \{x y z\}\]"
  puts $tcl_file "\$s2 set occupancy 0"
  puts $tcl_file "\$s1 set occupancy 1"
  puts $tcl_file "\$s2 writepdb initial_adjust.pdb"

  puts $tcl_file "mol delete all"
  puts $tcl_file "resetpsf"
  puts $tcl_file "mol load psf final_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_finmin/final_minimized\$cycle.coor"
  puts $tcl_file "set s1 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set s2 \[atomselect top \"all\"\]"
  puts $tcl_file "mol load pdb initial_target.pdb"
  puts $tcl_file "mol addfile starting_final.pdb_final_target.pdb_final_structure.dcd"
  puts $tcl_file "set s3 \[atomselect top \"name CA\"\]"
  puts $tcl_file "set trans_mat \[measure fit \$s1 \$s3\]"
  puts $tcl_file "\$s3 move \$trans_mat"
  puts $tcl_file "\$s1 set \{x y z\} \[\$s3 get \{x y z\}\]"
  puts $tcl_file "\$s2 set occupancy 0"
  puts $tcl_file "\$s1 set occupancy 1"
  puts $tcl_file "\$s2 writepdb final_adjust.pdb"

  puts $tcl_file "set sh_file \[open \"$output_prefix.sh\" w\]"
  puts $tcl_file "set sh_filename \"${output_prefix}.sh\""
  puts $tcl_file "puts \$sh_file \"\\\#\\\!\\\/bin\\\/bash\""
  puts $tcl_file "puts \$sh_file \"NAMD=\\\"\$namd2path \+p[expr ${num_cores}/2]\\\"\""
  puts $tcl_file "set namd_file \[open \[file join \"${output_prefix}_inipro\" \"pro.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ..\/initial_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ..\/initial_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname initial_process\$\{cycle\}\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  puts $tcl_file "puts \$namd_file \"parameters ../parameter_file.prm\""
  puts $tcl_file "puts \$namd_file \"set restartname res\""
  puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${output_prefix}_inimin\/initial_minimized\$\{cycle\}.coor\""
  puts $tcl_file "puts \$namd_file \"binvelocities ..\/${output_prefix}_inimin\/initial_minimized\$\{cycle\}.vel\""
  puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${output_prefix}_inimin\/initial_minimized\$\{cycle\}.xst\""
  puts $tcl_file "puts \$namd_file \"wrapWater on\""
  puts $tcl_file "puts \$namd_file \"wrapAll on\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"PME yes\""
  puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"TMD on\""
  puts $tcl_file "puts \$namd_file \"TMDk ${spring_k}\""
  puts $tcl_file "puts \$namd_file \"TMDOutputFreq 2500\""
  puts $tcl_file "puts \$namd_file \"TMDFile ..\/initial_adjust.pdb\""
  puts $tcl_file "puts \$namd_file \"TMDFirstStep 0\""
  puts $tcl_file "puts \$namd_file \"TMDLastStep [expr $tmd_len*500]\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies 2500\""
  puts $tcl_file "puts \$namd_file \"outputPressure 2500\""
  puts $tcl_file "puts \$namd_file \"restartfreq 2500\""
  puts $tcl_file "puts \$namd_file \"dcdfreq 2500\""
  puts $tcl_file "puts \$namd_file \"xstfreq 2500\""
  puts $tcl_file "puts \$namd_file \"run [expr $tmd_len*500]\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${output_prefix}_inipro\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD pro.conf > pro.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\""

  puts $tcl_file "set namd_file \[open \[file join \"${output_prefix}_finpro\" \"pro.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ../final_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ../final_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname final_process\$\{cycle\}\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  puts $tcl_file "puts \$namd_file \"parameters ../parameter_file.prm\""
  puts $tcl_file "puts \$namd_file \"set restartname res\""
  puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${output_prefix}_finmin\/final_minimized\$\{cycle\}.coor\""
  puts $tcl_file "puts \$namd_file \"binvelocities ..\/${output_prefix}_finmin\/final_minimized\$\{cycle\}.vel\""
  puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${output_prefix}_finmin\/final_minimized\$\{cycle\}.xst\""
  puts $tcl_file "puts \$namd_file \"wrapWater on\""
  puts $tcl_file "puts \$namd_file \"wrapAll on\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"PME yes\""
  puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"TMD on\""
  puts $tcl_file "puts \$namd_file \"TMDk ${spring_k}\""
  puts $tcl_file "puts \$namd_file \"TMDOutputFreq 2500\""
  puts $tcl_file "puts \$namd_file \"TMDFile ..\/final_adjust.pdb\""
  puts $tcl_file "puts \$namd_file \"TMDFirstStep 0\""
  puts $tcl_file "puts \$namd_file \"TMDLastStep [expr $tmd_len*500]\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies 2500\""
  puts $tcl_file "puts \$namd_file \"outputPressure 2500\""
  puts $tcl_file "puts \$namd_file \"restartfreq 2500\""
  puts $tcl_file "puts \$namd_file \"dcdfreq 2500\""
  puts $tcl_file "puts \$namd_file \"xstfreq 2500\""
  puts $tcl_file "puts \$namd_file \"run [expr $tmd_len*500]\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${output_prefix}_finpro\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD pro.conf > pro.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\""
  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"

  puts $tcl_file "set status \[catch \{exec prody catdcd initr.dcd \/${output_prefix}_inipro\/initial_process\$\{cycle\}.dcd -o initial_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv initial_trajectory.dcd initr.dcd\} output\]" 
  puts $tcl_file "set status \[catch \{exec prody catdcd fintr.dcd \/${output_prefix}_finpro\/final_process\$\{cycle\}.dcd -o final_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv final_trajectory.dcd fintr.dcd\} output\]" 
  
  ###### FINAL MINIMIZATION ########
  puts $tcl_file "set sh_file \[open \"$output_prefix.sh\" w\]"
  puts $tcl_file "set sh_filename \"${output_prefix}.sh\""
  puts $tcl_file "puts \$sh_file \"\\\#\\\!\\\/bin\\\/bash\""
  puts $tcl_file "puts \$sh_file \"NAMD=\\\"\$namd2path \+p[expr ${num_cores}/2]\\\"\""
  puts $tcl_file "set namd_file \[open \[file join \"${output_prefix}_inimin\" \"min.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ../initial_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ../initial_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname initial_minimized\[expr \$\{cycle\}+1\]\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  puts $tcl_file "puts \$namd_file \"parameters ../parameter_file.prm\""
  puts $tcl_file "puts \$namd_file \"set restartname res\""
  puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${output_prefix}_inipro\/initial_process\$\{cycle\}.coor\""
  puts $tcl_file "puts \$namd_file \"binvelocities ..\/${output_prefix}_inipro\/initial_process\$\{cycle\}.vel\""
  puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${output_prefix}_inipro\/initial_process\$\{cycle\}.xst\""
  puts $tcl_file "puts \$namd_file \"wrapWater on\""
  puts $tcl_file "puts \$namd_file \"wrapAll on\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"PME yes\""
  puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"outputPressure [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"restartfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"minimize [expr $min_length*500]\""
  puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${output_prefix}_inimin\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD min.conf > min.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\""

  puts $tcl_file "set namd_file \[open \[file join \"${output_prefix}_finmin\" \"min.conf\"\] w\]"
  puts $tcl_file "puts \$namd_file \"coordinates     ../final_ionized.pdb\""
  puts $tcl_file "puts \$namd_file \"structure       ../final_ionized.psf\""
  puts $tcl_file "puts \$namd_file \"set temperature $temperature\""
  puts $tcl_file "puts \$namd_file \"set outputname final_minimized\[expr \$\{cycle\}+1\]\""
  puts $tcl_file "puts \$namd_file \"set firsttimestep 0\""
  puts $tcl_file "puts \$namd_file \"paraTypeCharmm  on\""
  puts $tcl_file "puts \$namd_file \"parameters ../parameter_file.prm\""
  puts $tcl_file "puts \$namd_file \"set restartname res\""
  puts $tcl_file "puts \$namd_file \"bincoordinates ..\/${output_prefix}_finpro\/final_process\$\{cycle\}.coor\""
  puts $tcl_file "puts \$namd_file \"binvelocities ..\/${output_prefix}_finpro\/final_process\$\{cycle\}.vel\""
  puts $tcl_file "puts \$namd_file \"extendedSystem ..\/${output_prefix}_finpro\/final_process\$\{cycle\}.xst\""
  puts $tcl_file "puts \$namd_file \"wrapWater on\""
  puts $tcl_file "puts \$namd_file \"wrapAll on\""
  puts $tcl_file "puts \$namd_file \"exclude         scaled1-4\""
  puts $tcl_file "puts \$namd_file \"1-4scaling 1.0\""
  puts $tcl_file "puts \$namd_file \"cutoff 12.0\""
  puts $tcl_file "puts \$namd_file \"timestep        1.0\""
  puts $tcl_file "puts \$namd_file \"switching       on\""
  puts $tcl_file "puts \$namd_file \"switchdist      10.0\""
  puts $tcl_file "puts \$namd_file \"pairlistdist    13.5\""
  puts $tcl_file "puts \$namd_file \"rigidBonds none\""
  puts $tcl_file "puts \$namd_file \"nonbondedFreq 1\""
  puts $tcl_file "puts \$namd_file \"fullElectFrequency 1\""
  puts $tcl_file "puts \$namd_file \"stepspercycle 5\""
  puts $tcl_file "puts \$namd_file \"PME yes\""
  puts $tcl_file "puts \$namd_file \"PMEGridSpacing 1.0\"" 
  puts $tcl_file "puts \$namd_file \"langevin on\""
  puts $tcl_file "puts \$namd_file \"langevinDamping 5\""
  puts $tcl_file "puts \$namd_file \"langevinTemp \\\$temperature\""
  puts $tcl_file "puts \$namd_file \"langevinHydrogen on\""
  puts $tcl_file "puts \$namd_file \"outputname \\\$outputname\""
  puts $tcl_file "puts \$namd_file \"restartname \\\$restartname\""
  puts $tcl_file "puts \$namd_file \"outputEnergies [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"outputPressure [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"restartfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"dcdfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"xstfreq [expr $min_length*100]\""
  puts $tcl_file "puts \$namd_file \"minimize [expr $min_length*500]\""
  puts $tcl_file "puts \$namd_file \"reinitvels \\\$temperature\""
  puts $tcl_file "close \$namd_file"
  puts $tcl_file "puts \$sh_file \"cd ${output_prefix}_finmin\""
  puts $tcl_file "puts \$sh_file \"\\\$NAMD min.conf > min.log \&\""
  puts $tcl_file "puts \$sh_file \"cd ..\""
  puts $tcl_file "puts \$sh_file \"wait\""
  puts $tcl_file "close \$sh_file"
  puts $tcl_file "set status \[catch \{exec bash \$sh_filename\} output\]"

  puts $tcl_file "set status \[catch \{exec prody catdcd initr.dcd \/${output_prefix}_inimin\/initial_minimized\[expr \$\{cycle\}+1\].dcd -o initial_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv initial_trajectory.dcd initr.dcd\} output\]" 
  puts $tcl_file "set status \[catch \{exec prody catdcd fintr.dcd \/${output_prefix}_finmin\/final_minimized\[expr \$\{cycle\}+1\].dcd -o final_trajectory.dcd\} output\]"
  puts $tcl_file "set status \[catch \{exec mv final_trajectory.dcd fintr.dcd\} output\]" 
  
  puts $tcl_file "mol delete all" 
  puts $tcl_file "mol load psf initial_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_inimin/initial_minimized\[expr \$\{cycle\}+1\].coor" 
  puts $tcl_file "set sel1 \[atomselect top \"name CA\"\]" 
  puts $tcl_file "set sel1a \[atomselect top all\]"
  puts $tcl_file "mol load psf final_ionized.psf"
  puts $tcl_file "mol addfile ${output_prefix}_finmin/final_minimized\[expr \$\{cycle\}+1\].coor"  
  puts $tcl_file "set sel2 \[atomselect top \"name CA\"\]" 
  puts $tcl_file "set sel2a \[atomselect top all\]"
  puts $tcl_file "set trans_mat \[measure fit \$sel2 \$sel1\]"
  puts $tcl_file "\$sel2a move \$trans_mat"
  puts $tcl_file "set rmsd \[measure rmsd \$sel2 \$sel1\]"
  puts $tcl_file "set all_rmsd(\[expr \$\{cycle\}+1\]) rmsd"
  puts $tcl_file "puts \$rmsd"
  puts $tcl_file "if \{\(\$rmsd < 1.5)\|\|(\[expr \$all_rmsd\(\$\{cycle\}\) - \$all_rmsd\(\[expr \$\{cycle\}+1\]\)\]\ < 0.15 \)\} \{ break \}"
  puts $tcl_file "}"
  puts $tcl_file "set status \[catch \{exec mv initr.dcd initial_trajectory.dcd\} output\]" 
  puts $tcl_file "set status \[catch \{exec mv fintr.dcd final_trajectory.dcd\} output\]" 
  close $tcl_file
  file delete prop.pdb
  file delete pro.pdb
  file delete pro.psf
  file delete pro_wb.psf
  file delete pro_wb.pdb
  file delete pro_wb.log
  file delete fino.pdb
  file delete init.pdb
  file copy -force $topo_file $outputdir/topology_file.top
  file copy -force $para_file $outputdir/parameter_file.prm
  file copy -force $COMD_PATH/anmmc.py $outputdir/anmmc.py

  #source $tcl_file_name


  ::comd::Logview [file join "$outputdir" "$output_prefix.log"]

  tk_messageBox -type ok -title "Setup Complete" \
    -message "Setup of $output_prefix is complete. See $output_prefix.log file."


}

proc comd_tk {} {
  ::comd::comdgui
  #set ::druggability::which_mode [lindex $::druggability::titles [lsearch $::druggability::interfaces $::druggability::interface]]
  #::druggability::Switch_mode $::druggability::interface
  #return $::druggability::w
}



# proc drugui {args} {
#   global errorInfo errorCode
#   set oldcontext [psfcontext new]  ;# new context
#   set errflag [catch { eval drugui_core $args } errMsg]
#   set savedInfo $errorInfo
#   set savedCode $errorCode
#   psfcontext $oldcontext delete  ;# revert to old context
#   if $errflag { error $errMsg $savedInfo $savedCode }
# }


