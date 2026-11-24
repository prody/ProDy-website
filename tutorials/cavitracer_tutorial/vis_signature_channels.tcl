# VMD: signature channels and sequence conservation (CYP3A4/P450 ensemble)
# Run from the directory containing ensemble_out_channels/:
#     vmd -e vis_signature_channels.tcl
# The data come from calcSignatureChannels().

set result_dir "ensemble_out_channels"
set signature_file [file join $result_dir "signature_channel_pt20.pdb"]
set entropy_file   [file join $result_dir "signature_channel_entropy.csv"]

if {![file exists $signature_file]} {
    error "Missing $signature_file; run calcSignatureChannels first, or change result_dir."
}
if {![file exists $entropy_file]} {
    error "Missing $entropy_file; run the signature analysis first."
}

# Reference residue conservation is obtained from the entropy CSV.
# Do not use raw protein B-factors: residues lacking entropy retain the
# original crystallographic B-factor in the combined signature PDB.
array set entropy {}
set min_entropy 1.0e20
set max_entropy -1.0e20
set f [open $entropy_file r]
gets $f  ;# header: resnum,resname,entropy
while {[gets $f line] >= 0} {
    set fields [split $line ","]
    if {[llength $fields] < 3} { continue }
    set rn [string trim [lindex $fields 0]]
    set val [string trim [lindex $fields 2]]
    if {![string is integer -strict $rn] || ![string is double -strict $val]} { continue }
    set entropy($rn) [expr {double($val)}]
    if {$entropy($rn) < $min_entropy} { set min_entropy $entropy($rn) }
    if {$entropy($rn) > $max_entropy} { set max_entropy $entropy($rn) }
}
close $f

set mid [mol new $signature_file type pdb waitfor all]
mol delrep 0 $mid

set protein [atomselect $mid "protein"]
set conservation {}
set n_mapped 0
set n_missing 0
foreach rn [$protein get resid] {
    if {[info exists entropy($rn)]} {
        # Normalize the observed entropy range, then reverse it:
        # 1 (red) = greatest conservation, 0 (blue) = greatest variability.
        if {$max_entropy > $min_entropy} {
            set score [expr {($max_entropy - $entropy($rn)) / ($max_entropy - $min_entropy)}]
        } else {
            set score 1.0
        }
        lappend conservation $score
        incr n_mapped
    } else {
        # -1 marks residues absent from the sequence alignment.
        lappend conservation -1.0
        incr n_missing
    }
}
$protein set user $conservation
$protein delete

color scale method BWR

# Representation 0: protein NewCartoon, colored by relative conservation.
mol selection {protein and user >= 0}
mol representation NewCartoon
mol color User
mol material Opaque
mol addrep $mid
mol scaleminmax $mid 0 0.0 1.0

# Representation 1: residues with no usable Shannon entropy: grey.
mol selection {protein and user < 0}
mol representation NewCartoon
mol color ColorID 2
mol material Opaque
mol addrep $mid

# Representation 2: signature voxels as a smooth QuickSurf surface.
# FIL B-factors contain occupancy as percent (not Shannon entropy).
set voxels [atomselect $mid "resname FIL"]
if {[$voxels num] > 0} {
    set beta_values [$voxels get beta]
    set min_occ [lindex $beta_values 0]
    set max_occ [lindex $beta_values 0]
    foreach b $beta_values {
        if {$b < $min_occ} { set min_occ $b }
        if {$b > $max_occ} { set max_occ $b }
    }
    # Normalize to the actual occupancy range so the highest observed
    # occupancy is red even when no voxel reaches 100%.
    if {$max_occ <= $min_occ} { set max_occ [expr {$min_occ + 1.0}] }

    mol selection {resname FIL}
    mol representation QuickSurf 1.0 0.5 1.0 1.0
    mol color Beta
    mol material Transparent
    mol addrep $mid
    mol scaleminmax $mid 2 $min_occ $max_occ
    puts "Signature channels: FIL occupancy (blue = $min_occ%, red = $max_occ%)."
} else {
    puts "No FIL signature voxels found in $signature_file."
}
$voxels delete

axes location Off
display projection Orthographic
display resetview

puts "Protein conservation: BWR (blue = variable, red = conserved)."
puts "Mapped protein atoms: $n_mapped; atoms without entropy: $n_missing (grey)."
puts "Sequence entropy range: $min_entropy to $max_entropy."
