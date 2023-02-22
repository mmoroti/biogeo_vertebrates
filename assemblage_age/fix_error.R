#######################################################
# This script is a rough-and-ready attempt at 
# checking a list of basic issues that people often 
# have with their input tree files and geography files.
# 
# In the case of trees, some of the issues have 
# standard fixes.
#
# Mostly I am posting this to save myself time, when
# diagnosing peoples' issues, but feel free to use 
# it and save me the time!
#######################################################

# Packages/updates to source
library(ape)
library(BioGeoBEARS)
sourceall("/drives/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
# (you could alternatively use the source() commands in 
#  the example script, at:
#  http://phylo.wikidot.com/biogeobears#script
# )

# Set your working directory
wd = "/drives/GDrive/__GDrive_projects/2016-02-09_Robin_Beck/_v4/_run6_BDSS_M0_7areas/"
setwd(wd)

# Set your tree file and geography file names
#trfn = "tree5.newick"
trfn = "assemblage_age/mammals_biogeo_harm.new"
geogfn = "assemblage_age/geo_area_mammals_harm.data"

#######################################################
# Load the tree
# 
# Your tree file should be Newick format. If, for 
# some reason, you've got a NEXUS file, this will attempt to 
# convert to Newick.
#######################################################

# Tree file
try_result = try(read.tree(trfn))
nexTF = grepl(pattern="missing value where TRUE/FALSE needed", x=try_result)
if ( (length(nexTF) == 1) && (nexTF == TRUE) )
{
  cat("\n\nERROR: Reading the tree file as NEXUS resulted in this error:\n\n")
  print(try_result)
  cat("\n...probably this means it is NEXUS rather than Newick. We will attempt to\nread it as NEXUS next, and save to Newick.")
  try_result = try(read.nexus(trfn))
  if (class(try_result) == "phylo")
  {
    cat("\n\nFile was NEXUS. Converting to Newick.\n\n")
    new_trfn = paste0(trfn, ".newick")
    cat(trfn, " --> ", new_trfn)
    write.tree(tr, file=new_trfn)
    trfn = new_trfn
    tr = try_result
  } # END if (class(try_result) == "phylo")
  if (class(try_result) == "try-error")
  {
    stoptxt = "ERROR: Reading the tree failed for both read.tree() and read.nexus. Look at your file and verify the format!"
    cat("\n\n")
    cat(stoptxt)
    cat("\n\n")
    stop()
  } else {
    tr = try_result
  }# END if (class(try_result) = "try-error")
} else {
  # Presumably it's Newick
  tr = try_result
} # END if (nexTF == TRUE)

plot(tr)
title("Newick tree")
axisPhylo()

# Is the tree dichotomous?  And rooted?
if (is.binary.tree(tr) == FALSE)
{
  stop("Stopping, because tree is not binary. (In other words, you have polytomies. BioGeoBEARS cannot handle polytomies.)")
}
if (is.rooted(tr) == FALSE)
{
  stop("Stopping, because tree is not rooted. (BioGeoBEARS requires rooted trees.)")
}

#######################################################
# Check for negative and 0-length branches
# Note: in BioGeoBEARS, tips with branches very close to length 0 (but not exactly 0)
#       are treated as direct ancestors.  But, negative branchlengths, and internal
#       branchlengths=0, all have to be fixed.
#######################################################
# Minimum branchlength for any new internal branches
min_brlen = 0.01

# Are there 0-length branches?
min(tr$edge.length)

# Remove branchlengths less than 0
if (min(tr$edge.length) < 0)
{
  cat("\n\n")
  cat("Negative branchlengths detected at these nodes:\n\n")
  trtable = prt(tr, printflag=FALSE)
  negative_brlens_TF = trtable$edge.length < 0
  negative_brlens_TF[is.na(negative_brlens_TF)] = FALSE
  negative_BL_rows = trtable[negative_brlens_TF, ]
  print(negative_BL_rows)
  cat("\n")
  
  cat("Correcting, using impose_min_brlen()")
  tr = impose_min_brlen(phy=tr, min_brlen=min_brlen, leave_BL0_terminals=TRUE)
  
  if (grepl(pattern="\\.newick", x=trfn) == TRUE)
  {
    new_trfn = gsub(pattern="\\.newick", replacement="_noNeg.newick", x=trfn)
  } else {
    new_trfn = paste0(trfn, "_noNeg.newick")
  }
  trfn = new_trfn
  
  write.tree(tr, file=new_trfn)
  cat("\n\nFixed negative branchlengths, and saved to ", new_trfn, "\n", sep="")
} # END if (min(tr$edge.length) < 0)
min(tr$edge.length)
sum(tr$edge.length == min(tr$edge.length))

#######################################################
# *Internal* branches of 0 length?
#######################################################
trtable = prt(tr, printflag=FALSE)
internal_TF = trtable$node.type == "internal"
edges_BL0_TF = trtable$edge.length == 0
sum_TFs = (internal_TF + edges_BL0_TF)
sum_TFs[is.na(sum_TFs)] = 0
internal_BL0_TF = (sum_TFs == 2)
sum(internal_BL0_TF)

# Edit the branchlengths, if needed
if (sum(internal_BL0_TF) > 0)
{
  internal_BL0_rows = trtable[internal_BL0_TF, ]
  internal_BL0_rows
  nodes_to_change = internal_BL0_rows$node
  edges_to_change = internal_BL0_rows$parent_br
  edges_to_change
  
  cat("\n\n")
  cat("Internal branches of length 0 detected at these nodes:\n\n")
  print(internal_BL0_rows)
  cat("\n")
  cat("Changing to length min_brlen=", min_brlen, "...")
  tr$edge.length[edges_to_change] = min_brlen
  
  # New filename
  if (grepl(pattern="\\.newick", x=trfn) == TRUE)
  {
    new_trfn = gsub(pattern="\\.newick", replacement="_minBL.newick", x=trfn)
  } else {
    new_trfn = paste0(trfn, "_minBL.newick")
  }
  trfn = new_trfn
  
  write.tree(tr, file=new_trfn)
  cat("...saved to ", new_trfn, "\n", sep="")
} # END if (sum(internal_BL0_TF) > 0)

##################################################
# *External* branches of 0 length?
##################################################
external_TF = trtable$node.type == "tip"
edges_BL0_TF = trtable$edge.length == 0
sum_TFs = (external_TF + edges_BL0_TF)
sum_TFs[is.na(sum_TFs)] = 0
external_BL0_TF = (sum_TFs == 2)
sum(external_BL0_TF)

tip_brlen = 0.0000001
if (sum(external_BL0_TF) > 0)
{
  cat("\n\n")
  cat("Tip branches of length zero detected. Changing to length ", tip_brlen, sep="")
  cat("\n\n")
  edgenums_to_change = trtable$parent_br[external_BL0_TF]
  edgenums_to_change
  
  # Change them
  tr$edge.length[edgenums_to_change] = tip_brlen
  
  # New filename
  if (grepl(pattern="\\.newick", x=trfn) == TRUE)
  {
    new_trfn = gsub(pattern="\\.newick", replacement="_tipsNo0.newick", x=trfn)
  } else {
    new_trfn = paste0(trfn, "_tipsNo0.newick")
  }
  trfn = new_trfn
  
  write.tree(tr, file=new_trfn)
  cat("...saved to ", new_trfn, "\n", sep="")
} # END if (sum(external_BL0_TF) > 0)

#######################################################
# Try loading the geography file
#######################################################
# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

#######################################################
# Final check
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Specify tree file and geography file
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
