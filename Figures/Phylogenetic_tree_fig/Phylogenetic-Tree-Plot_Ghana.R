# Load libraries
library(ggtree)
library(ggplot2)
library(treeio)
library(dplyr)
library(tidytree)
library(jsonlite)
library(RColorBrewer)
library(cowplot)
library(rlang)
library(ape)
library(phytools)

# Sets the script's working directory to its current location. By
# doing this, file paths can be written in reference to the scripts location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#' Function:
#' @name parseNSCladesJSON
#' @description Transforms the Clade JSON structure produced from Nextstrain into a
#'     dataframe that is linked to the actual node name from the tree data. The JSON file
#'     uses Node text labels rather than NodeIDs, which can be used to query the tree object itself.
#'     Thus, creating a dataframe that links a more useful identifier to
#'     the clade information is helpful for subsetting and collapsing clades on the tree. Additionally,
#'     this new dataframe format makes it useful for identifying which nodes are the start of a given clade.
#'     In the Clade JSON file from Nextstrain, nodes that represent the start of a new clade contain a
#'     "clade_annotation" attribute. Thus, in a dataframe format, we can easily identify nodes that
#'     are the top of a clade by finding those that have a value for this column. In a new nextstrain update,
#'     the format of the cladeJSON file has been modified such that the "clade_annotation" designation does
#'     not exist, but instead a "branch" designation is assigned to those noes that are indeed the branches.
#'     This function has been modified to import the new cladeJSON format, yet map it in such a way that if 
#'     the node is listed as a branch, then the lineage information is filled in for "clade_annotation" so that 
#'     the downstream code would still work.
#'
#' @param cladeJSON A JSON file produced by Nextstrain containing information regarding what clade each
#'      node (or sample) falls under on the phylogenetic tree.
#'      
#' @param treeData A Dataframe representation of a phylogenetic tree imported through
#'      the treeIO package.
#' 
#' @return a dataframe containing columns Node, label, clade_membership, outbreak, lineage,
#'       clade_annotation
#'       
parseNSCladesJSON <- function(cladeJSON, treeData) {
  # Create an empty dataframe to store the parsed clade data
  cladeDF <- data.frame(node = character(), label = character(), 
                        clade_membership = character(), outbreak = character(), 
                        lineage = character(), clade_annotation = character(), 
                        stringsAsFactors = FALSE)
  
  # Loops over each node in the clade JSON file
  if ("nodes" %in% names(cladeJSON)) {
    cladeNodes <- names(cladeJSON$nodes)
    
    for (i in cladeNodes) {
      node <- ifelse(i %in% treeData$label, treeData$node[treeData$label == i], NA)
      cm <- cladeJSON$nodes[[i]]$clade_membership
      ob <- cladeJSON$nodes[[i]]$outbreak
      ln <- cladeJSON$nodes[[i]]$lineage
      
      # Default clade annotation to NA
      ca <- NA
      
      # If node is listed in branches and has a clade label, update clade annotation
      if ("branches" %in% names(cladeJSON) && i %in% names(cladeJSON$branches)) {
        if ("labels" %in% names(cladeJSON$branches[[i]]) && "clade" %in% names(cladeJSON$branches[[i]]$labels)) {
          ca <- cladeJSON$branches[[i]]$labels$clade
        }
      }
      
      # Append to the dataframe
      cladeDF <- rbind(cladeDF, data.frame(node = node, label = i, 
                                           clade_membership = cm, outbreak = ob, 
                                           lineage = ln, clade_annotation = ca))
    }
  }
  
  return(cladeDF)
}


#' Function:
#' @name getParentLineageNodes
#' @description Finds the Node IDs for the clades that a set of Node IDs 
#'      falls under and creates a new dataframe containing the combined data.
#'      This is different from simply grabbing the parent node. The direct
#'      parent node is simply the node above a given node. However, this node may
#'      not be a node defining a given clade.
#'
#' @param nodes A list of nodes.
#'      
#' @param treeData A Dataframe representation of a phylogenetic tree imported through
#'      the treeIO package.
#' 
#' @param lineageNodes A Dataframe containing metadata for the Nodes that represent 
#'      the top each clade.
#'
#'
#' @return a dataframe containing columns node, lineage, lineageParentNodeName, lineageParenteNodeNum
#'       
getParentLineageNodes <- function(nodes, treeData, lineageNodes) {
  
  # Creates an empty dataframe to store final data.
  # The dataframe will consist of the following columns:
  #   - node: The Node ID
  #   - lineage: the lineage of the node.
  #   - lineageParentNodeName: the name of the lineage's top node.
  #   - lineageParentNodeNum: the nodeID of the lineage's top node.
  returnDF <- data.frame(matrix(ncol=4, nrow=0))
  colnames(returnDF) <- c("node", "lineage", "lineageParentNodeName", "lineageParentNodeNum")
  
  # Loop over the list of nodes provided.
  for (n in nodes) {
    print(paste0("Checking Node: ", n))
    
    # To find the lineage parent for a given node, we can use a recursive style method,
    # where we find the continue identifying individual parent nodes until a 
    # lineage parent node is identified. Once this is found, we can end the recursion.
    
    # Sets a boolean value to denote once a lineage parent node has been found.
    lineageNodeFound <- FALSE
    
    # Sets the current node being searched to the current node in the list.
    nodeToSearch <- n
    
    # While a lineage parent has not been found.
    while (!lineageNodeFound) {
      # Grab the parent node and node label for the current
      # node being searched.
      parentNode <- treeDF$parent[treeDF$node == nodeToSearch]
      parentLabel <- treeDF$label[treeDF$node == parentNode]
      
      #print(parentNode)
      #print(parentLabel)
      
      # Check whether the parent node is one of the lineage
      # parents.
      if (parentNode %in% as.list(lineageNodes$node)) {
        
        # If yes, then set the boolean value to true to end
        # the recursion.
        lineageNodeFound <- TRUE
        
        # Grab the lineage value of the parent node
        lin <- lineageNodes$lineage[lineageNodes$node == parentNode]
        print(paste0("This node is was found to be a lineage node of lineage: ", lin))
        print(lin)
        
        # Append a row to the return database to include the base node, the lineage, the
        # lineage parent node label, and lineage parent node ID.
        row <- data.frame(node=n, lineage=lin, lineageParentName=parentLabel, lineageParentNode=parentNode)
        returnDF <- rbind(returnDF, row)
      }
      else {
        # If not, set the node to be searched equal to the current parent node ID.
        nodeToSearch <- parentNode 
      }
    }
  }
  
  # Once all of the nodes have been processed, return the dataframe.
  return(returnDF)
}


#' Function:
#' @name getNodesNecessaryToDisplay
#' @description Finds the minimum node IDs necessary to display a set of nodes
#'
#' @param nodes A list of nodes.
#'      
#' @param treeData A Dataframe representation of a phylogenetic tree imported through
#'      the treeIO package.
#'
#' @return a list of unique nodes
#'  
getNodesNecessaryToDisplay <- function(nodes, treeData) {
  
  # Creates an empty list to store the nodes to be
  # returned.
  toReturn <- c()
  
  # Loop over the list of nodes.
  for (n in nodes) {
    
    # Since we are trying to determine which nodes are necessary to
    #'     display tha given node. We can recursively traverse up
    # the tree structure until the root is reached. By doing so, 
    # we will have found a list of nodes that go from the root to
    # the given leaf node.
    
    # Create a boolean value to act as the stop condition
    # of the recursion when the root is reached.
    reachedRoot = FALSE
    
    # Set the current node equal to the one being checked.
    currNode <- n
    
    # Loop until the root node has been reached.
    while (!reachedRoot) {
      # If we have not reached the root node yet,
      # add the current node to the list of nodes
      toReturn <- append(toReturn, currNode)
      
      # Grab the parent for the current node.
      parent <- treeDF$parent[treeDF$node == currNode]
      
      # The root node has an interesting property in that the
      # the "parent node ID is the same as its individual ID. Thus,
      # we can use this do determine whether we have reached the root node.
      if (parent == currNode) {
        # If yes, set hte boolean value to true.
        reachedRoot = TRUE
      }
      else {
        # If not, set the current node equal to the
        # parent node so that it will be processed during
        # the next recursive iteration.
        currNode = parent
      }
    }
  }
  
  # Once all of the input nodes have been processed, remove
  # duplicate nodes from the master list of nodes necessary
  # to display each input node and return this list. 
  return(unique(toReturn))
}


#' Function:
#' @name collapseMultipleNodes
#' @description Collapses a set of nodes on a tree plot. ggtree only allows
#'     for one node to be collapsed at a time. Thus, this command uses
#'     a for loop to iterate over the individual nodes and collapse
#'     them one by one, saving the tree each time. The collapsed nodes
#'     are represented using a black diamon shape on the tree. 
#'
#' @param plot A ggtree plot object to be modified.
#'      
#' @param nodes A list of nodes to collapse in the plot.
#'
#' @return A modified tree where all of the provided nodes have been collapsed.
#'  
collapseMultipleNodes <- function(plot, nodes) {
  
  # Identify the root node
  root_node <- plot$data$node[plot$data$node == plot$data$parent]
  print(root_node)
  
  # Track collapsed nodes
  collapsed_nodes <- integer(0)
  print(collapsed_nodes)
  
  # Loops over the nodes in the provided list.
  for (n in nodes) {
    print(paste("Processing node:", n))  # Debugging information
    
    # Check if the node is already in the list of collapsed nodes
    if (n %in% collapsed_nodes) {
      print(paste("Node", n, "is already collapsed"))
      next
    }
    
    # Collapses the current node on the tree using the provided collapse function, and then
    # adds a black diamond shape to the tree at the position that the node was collapsed.
    plot <- collapse(plot, node = n) + geom_point2(aes(subset = (node == n)), shape = 23, size = 2, fill = 'black')
    
    # Add the node to the list of collapsed nodes
    collapsed_nodes <- c(collapsed_nodes, n)
    print(collapsed_nodes)
  }
  
  # Adds the shapes again (Need to revisit this code to check whether necessary).
  plot <- plot + geom_point2(aes(subset = (node %in% nodes)), shape = 23, size = 2, fill = 'black')
  
  # Returns the modified plot.
  return(plot)
}


#' Function:
#' @name addCladeLabels
#' @description Adds a list of clade labels to a tree plot. The geom_cladelab
#'      function provided by ggtree only allows for the addition of one label
#'      at a time. Thus, this function uses iteration to add each label one
#'      by one and returns the final plot.
#'
#' @param plot A ggtree plot object to be modified.
#'      
#' @param nodeDF A Dataframe containing information about a set of nodes and
#'      and what label to give them.
#'
#' @return A modified tree where all of the provided nodes have been labeled.
#'  
addCladeLabels <- function(plot, nodeDF) {
  
  # Loop over each node in the dataframe.
  for (n in nodeDF$nodes) {
    # Add a geom_cladelab to the current node, where the label
    # is the corresponding 'label' attribute in the dataframe, to
    # the tree.
    plot <- plot + geom_cladelab(node = n, label=nodeDF$label[nodeDF$node == n])
  }
  
  # Return the finalized tree.
  return(plot)
}

# Creates a list of pre-defined sample names for the Ghana Mpox samples.
GHSampleNames <- c("22-007", "22-010", "22-030", "22-039", 
                   "22-046", "22-106", "22-173", "22-177", 
                   "22-190", "22-195", "22-240", "22-244", "22-328", 
                   "22-345", "22-401", "22-467", "22-496A")


# Reads in the tree in Nexus format from a file (this can be downloaded from Auspice when viewing the tree produced
# by nextstrain).
tree <- read.nexus("../../../Analysis/Ghana_2024/MPXV/20241211_Ghana_MPXV_final_analysis/Nextstrain_output_r_analysis/nextstrain_Ghana_mpxv_all-cladeII_17_new_samples-final-7_tree_midpoint_rooted.nexus")

# Converts the tree object into a dataframe format.
treeDF <- as.data.frame(as_tibble(tree))

# Reads in the clade information produced by nextstrain in JSON format using the jsonlite package.
clades <- jsonlite::fromJSON("../../../Analysis/Ghana_2024/MPXV/20241211_Ghana_MPXV_final_analysis/Nextstrain_output_r_analysis/clades.json", flatten = TRUE)

# Combines the tree data and clade information into a singe dataframe.
# This will be used to subset and filter the tree based on clade
# and lineage assignments.
cladeDF <- parseNSCladesJSON(clades, treeDF)

# Subset the dataframe to obtain a dataframe containing only data
# for the nodes that 'define' each clade/lineage. This can be identified by
# choosing only nodes who have a value in the "clade_annotation" attribute. This
# is how nextstrain formats their JSON files. Secondly, for the GH Samples, if 
# the lineage field is blank, fill it from clade_annotation so that there are no
# blanks
lineageNodeDF <- cladeDF %>% 
  filter(!is.na(clade_annotation)) %>%
  mutate(lineage = ifelse(lineage == "", clade_annotation, lineage))
print(head(lineageNodeDF))
# Print unique values of the column "node" in lineageNodeDF
#unique_nodes <- unique(lineageNodeDF$node)
#print(unique_nodes)

# Print unique values of the column "node" in lineageNodeDF in a more readable format
#cat("Unique node values in lineageNodeDF:\n", paste(unique_nodes, collapse = ", "))


# Create a dataframe containing the tree data (node, lineage, etc) for only
# the Ghana samples by joining the Ghana list to the
# tree dataframe as a whole.
GHSampleTreeData <- data.frame(Sample=GHSampleNames)
GHSampleTreeData$node <- treeDF$node[match(GHSampleTreeData$Sample, treeDF$label)]


#getParentLineageNodes(GHSampleTreeData$node, treeDF, lineageNodeDF)

# Add in extra columns to this dataframe containing information about parent
# nodes of the lineage that the GH samples fall under. This can be
# done by joining the existing dataframe with the result of the getParentLineageNodes
# function on the "node" attribute.
GHSampleTreeData <- merge(GHSampleTreeData, getParentLineageNodes(GHSampleTreeData$node, treeDF, lineageNodeDF), by="node")

# Filter to keep only unique rows with non-NA values in the lineage column
GHSampleTreeData <- GHSampleTreeData %>%
  filter(!is.na(lineage)) %>%
  distinct()

print(head(GHSampleTreeData))

# Get a list of the lineage parent nodes that have no GH samples
# classified under them. These will be the nodes to collapse (as they do not contain any
# GH samples and will clutter the tree). This is accomplished by finding the difference in the
# set of all lineage parent nodes and the unique list of lineage parent nodes from the 
# GH sample data.
nonGHLineageParents <- setdiff(unique(lineageNodeDF$node), unique(GHSampleTreeData$lineageParentNode))
print(nonGHLineageParents)


# Grab a list of nodes that are necessary to display the GH Sample lineage parent nodes.
neededToDisplay <- getNodesNecessaryToDisplay(GHSampleTreeData$lineageParentNode, treeDF)
print(neededToDisplay)

#nodesToAddToDisplay <- c(1335)
#finalNeededToDisplay <- union(neededToDisplay, nodesToAddToDisplay)
#print(finalNeededToDisplay)

# Uses the setdiff function to find all of the nodes we can collapse without losing any
# that are required for the GH samples, excluding the nodes you want to preserve.
nodesToCollapse <- setdiff(treeDF$node[match(unique(lineageNodeDF$node), treeDF$node)], neededToDisplay)
print(nodesToCollapse)

#nodesToRemoveFromNodesToCollapse <- c(1226)
#finalNodesToCollapse <- setdiff(nodesToCollapse, nodesToRemoveFromNodesToCollapse)
#print(finalNodesToCollapse)

# Ensure treeDF has a column 'node' that contains the node IDs
sortedNodesToBeCollapsed <- treeDF %>%
  filter(node %in% nodesToCollapse) %>%
  arrange(match(node, nodesToCollapse)) %>%
  pull(node)

# Verify the order of nodes
print(sortedNodesToBeCollapsed)

# Find common nodes
common_nodes <- intersect(neededToDisplay, sortedNodesToBeCollapsed)
print(common_nodes)


# A test command to grab all of the nodes needed to collapse the B1 sublineages.
# This lineage might need to be changed based on the nextstrain run and how many
# samples are in the NS build. Every instance of this '1335' would need to be 
# updated for the tree to be generated with no errors. This is the position 
# of the B.1 clade (1335 in this dataset) and the A.1.1 clade (1334) its parent.

# Define the nodes you want to preserve
#nodesToPreserve <- lineageNodeDF$node[grepl("^A|clade|^B\\.1$|^B\\.1\\.2$", lineageNodeDF$lineage)]
#print(nodesToPreserve)

# Collapse all B.1 sublineages except the ones we want to preserve
#collapseAllB1Sublins <- setdiff(lineageNodeDF$node, finalNeededToDisplay)
#print(collapseAllB1Sublins)

# Grabs a dataframe of all of the GH sample points on the tree with
# their corresponding label. This will be used to label the GH samples
# on the tree only.
GHSamplePoints <- GHSampleTreeData %>% reframe(node, label=Sample)
print(head(GHSamplePoints))

# Grab a dataframe containing all of the clade nodes and their
# corresponding labels to be placed on the tree. The B1 label is also
# removed to simplify visualization.
GHSampleClades <- GHSampleTreeData %>% distinct(node=lineageParentNode, label=lineage)

# Assuming GHSampleClades is previously defined and contains the 'node' column
#GHSampleCladesMinusB1 <- GHSampleClades %>%
  #filter(node != 1335)
  

# Print the resulting data frame to verify
print(GHSampleClades)
#print(GHSampleCladesMinusB1)

#New code for adding isTip column to treeDF for plotting purposes:
# Assuming treeDF has columns 'node' and 'parent'
# Mark nodes as 'isParent' if they appear in the 'parent' column of any row
treeDF <- treeDF %>%
  mutate(isParent = node %in% parent)

# Nodes that are not listed as 'parent' anywhere are tips
treeDF <- treeDF %>%
  mutate(isTip = !isParent)

# Now remove the temporary 'isParent' column
treeDF <- select(treeDF, -isParent)

# To check the structure and make sure 'isTip' is created correctly
print(head(treeDF))
print(treeDF %>% filter(node %in% GHSamplePoints$node))
#print(treeDF %>% filter(node == root_node))

# Ensure `collapsed` column is correctly defined in treeDF
if (!"collapsed" %in% colnames(treeDF)) {
  treeDF$collapsed <- FALSE
}

# Ensure the 'x' and 'y' columns are initialized
if (!("x" %in% colnames(treeDF))) {
  treeDF$x <- NA
}
if (!("y" %in% colnames(treeDF))) {
  treeDF$y <- NA
}

# This code will add an 'isTip' column where TRUE indicates a tip node
# Now you can use this to filter tip nodes for plotting or other analyses
filtered_nodes <- treeDF %>%
  filter(isTip & node %in% GHSamplePoints$node)

# Print summary of the filtered nodes
print(summary(filtered_nodes))

# Ensure `isTip` column is correctly defined
if (!"isTip" %in% colnames(treeDF)) {
  treeDF <- treeDF %>%
    mutate(isTip = !node %in% parent)
}

# Ensure `collapsed` column is correctly defined in treeDF
if (!"collapsed" %in% colnames(treeDF)) {
  treeDF$collapsed <- NA
}

# Check if the `isTip` column and `GHSamplePoints` nodes are correctly identified
print(treeDF %>% filter(node %in% GHSamplePoints$node))
# Check if any subset used in plotting could be empty
print(nrow(treeDF %>% filter(isTip, node %in% GHSamplePoints$node)))

# Get pairwise distances between all nodes
distances <- cophenetic(tree)

# Calculate the total distance to all tips for each internal node
total_distances <- colSums(distances)

# Find the midpoint as the node with the smallest difference in distances
midpoint_node <- names(which.min(abs(total_distances - mean(total_distances))))

# Print the midpoint node
print(paste("The midpoint node is:", midpoint_node))

# Optionally, plot the tree and highlight the midpoint
#plot(tree, main = "Tree with Midpoint Highlighted")
#nodelabels(node = as.numeric(midpoint_node), pch = 19, col = "red")

print(nrow(treeDF %>% filter(!isTip, node == 169))) # midpoint root node from function above
#print(nrow(treeDF %>% filter(!isTip, node == 225))) # clade I node, not needed in this build as the tree is midpoint rooted
print(nrow(treeDF %>% filter(isTip, node == 696))) # lineage B.1.2 node
print(nrow(treeDF %>% filter(!isTip, node %in% GHSampleClades$node)))
# Check if node root exists and its status
print(treeDF %>% filter(node == 169))
# Check if node clade I node exists and its status
#print(treeDF %>% filter(node == 225)) # clade I node, not needed in this build as the tree is midpoint rooted
# Check if node lineage B.1.2. exists and its status
print(treeDF %>% filter(node == 696))



### Tree Plots

####################
## This first plot creates a tree with the following characteristics:
# Creates a general tree visualization with the root is labeled and marked with
# an open green circle. Additionally, the GH samples are marked with a red circle on the tree.

rawTreePlot <- ggtree(tree) + coord_cartesian(clip="off") +
  geom_tippoint(data=td_filter(isTip & node %in% GHSamplePoints$node), size=2, shape=21, color="black", fill="#228B22") +
  geom_tippoint(data=td_filter(isTip & node == 169), size=2, shape=21, color="black", fill="red") +
  geom_tiplab(data = td_filter(isTip & node == 169), size = 2.5, nudge_x = 0.00006) + 
  geom_treescale(x=0, y=45, offset=1)

# Print initial plot
print(rawTreePlot)


#
# Creates a general tree visualization with:
# - A red circle for the midpoint root at node == 169 (TRM291)
# - Green circles for the Ghana samples
# - A legend overlayed on the tree in the upper right corner

rawTreePlot2 <- ggtree(tree) + 
  coord_cartesian(clip = "off") +
  
  # Ghana samples (green circles with black outline)
  geom_tippoint(
    data = td_filter(isTip & node %in% GHSamplePoints$node), 
    aes(fill = "Ghana samples"),  # Use `fill` for legend mapping
    size = 2, shape = 21, color = "black"  # Outline color is explicitly black
  ) +
  
  # Midpoint root (red circle with black outline)
  geom_tippoint(
    data = td_filter(isTip & node == 169), 
    aes(fill = "Midpoint root (TRM291)"),  # Use `fill` for legend mapping
    size = 2, shape = 21, color = "black"  # Outline color is explicitly black
  ) +
  
  # Tree scale
  geom_treescale(x = 0, y = 45, offset = 1) +
  
  # Add a legend for fill colors
  scale_fill_manual(
    name = "Legend", 
    values = c(
      "Ghana samples" = "#228B22",  # Forest green
      "Midpoint root (TRM291)" = "red"
    ),
    labels = c(
      "Midpoint root (TRM291)" = "Midpoint root (TRM291)",
      "Ghana samples" = "Ghana samples"
    )
  ) +
  
  # Customize legend placement to overlay on the tree
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 4, 1, 1), "cm"),
    legend.position = c(0.7, 0.6),  # Upper right corner (x, y coordinates)
    legend.background = element_rect(fill = "white", color = "black", size = 0.5),  # Background for the legend
    legend.key = element_blank()  # Remove boxes around legend symbols
  ) +
  
  # Reduce the legend size
  guides(
    fill = guide_legend(override.aes = list(size = 3, shape = 21, color = "black"))  # Ensure legend uses black outline
  )

# Print the updated plot
print(rawTreePlot2)

# Remove the tree scale from rawTreePlot by recreating it without the tree scale. This is to avoid the treescale be duplicated in the 
# collapsed tree. 
rawTreePlot_no_scale <- ggtree(tree) + coord_cartesian(clip = "off") +
  geom_tippoint(data = td_filter(isTip & node %in% GHSampleTreeData$node), size = 2, shape = 21, color = "black", fill = "#228B22") +
  geom_tippoint(data = td_filter(isTip & node == 169), size = 2, shape = 21, color = "black", fill = "red") +
  geom_tiplab(data = td_filter(isTip & node == 169), size = 2.5, nudge_x = 0.00006)

print(rawTreePlot_no_scale)
# Collapse nodes and add the tree scale at the new position
treePlotCollapsed <- collapseMultipleNodes(rawTreePlot_no_scale, node = sortedNodesToBeCollapsed) +
  geom_treescale(x = 0, y = 40, offset = 1)

# Print final plot data and plot
print(treePlotCollapsed$data)
print(treePlotCollapsed)

# Adds a special label for the B.1 clade, as including this as a vertical
# line covers over some of the other clade labels. Then, clade labels for the other
# clades are added individually to customize position and spacing.
finalTreePlot <- treePlotCollapsed %<+% cladeDF + 
  # Label for lineage B.1 
  geom_nodelab(
    data=td_filter(!isTip & node == 647), 
    aes(label = "B.1"), 
    geom="text", 
    nudge_x = -0.000035, 
    nudge_y = 1,
    hjust= 0, 
    vjust= 0
  ) +
  # Label for lineage A.1.1
  geom_nodelab(
    data=td_filter(!isTip & node == 646), 
    aes(label=lineage), 
    geom="text", 
    nudge_x = -0.000055, 
    nudge_y = 1,
    hjust= 0, 
    vjust= 0
  ) +
  # Label for clade IIa
  geom_nodelab(
    data = td_filter(!isTip & node == 502), 
    aes(label = "Clade IIa"), 
    geom = "text", 
    nudge_x = 0.000945, 
    nudge_y = -14,
    hjust = 0, 
    vjust = 0
  ) +
  # Label for lineage A
  geom_nodelab(
    data = td_filter(!isTip & node == 528), 
    aes(label = "A"), 
    geom = "text", 
    nudge_x = -0.0000175, 
    nudge_y = 1,
    hjust = 0, 
    vjust = 0
  ) +
  # Label for lineage A.1
  geom_nodelab(
    data = td_filter(!isTip & node == 636), 
    aes(label = "A.1"), 
    geom = "text", 
    nudge_x = -0.000035, 
    nudge_y = 1,
    hjust = 0, 
    vjust = 0
  ) +
  # Label for clade IIb
  geom_nodelab(
    data = td_filter(!isTip & node == 547), 
    aes(label = "Clade IIb"), 
    geom = "text", 
    nudge_x = -0.000285, 
    nudge_y = -25,
    hjust = 0, 
    vjust = 0
  ) +
  # Label for clade I
  #geom_nodelab(
    #data = td_filter(!isTip & node == 225), 
    #aes(label = "Clade I"), 
    #geom = "text", 
    #nudge_x = -0.000065, 
    #nudge_y = -18,
    #hjust = 1, 
    #vjust = -2.75
  #) +
  # Label for clade A.2
  geom_nodelab(
    data = td_filter(!isTip & node == 561), 
    aes(label = "A.2"), 
    geom = "text", 
    nudge_x = -0.000055, 
    nudge_y = 1,
    hjust = -0.5, 
    vjust = 0
  ) +
  # Label for lineage A.2.1
  geom_nodelab(
    data = td_filter(!isTip & node == 576), 
    aes(label = "A.2.1"), 
    geom = "text", 
    nudge_x = 0.000205, # Increase this value slightly to move the label to the right
    nudge_y = 0,    # Adjust this value if you need to move the label up or down
    hjust = 0, 
    vjust = -1.35
  ) +
  geom_segment(
    data = td_filter(!isTip & node == 576), 
    aes(x = x + 0.0002, xend = x + 0.0002, y = y - 2, yend = y + 15), # Adjust y and yend to move the bar upwards and center it
    color = "black"
  ) +
  # Label for clade A.2.2, not needed b/c it is collapsed
  #geom_nodelab(
    #data = td_filter(!isTip & node == 563), 
    #aes(label = "A.2.2"), 
    #geom = "text", 
    #nudge_x = 0.000195, 
    #nudge_y = 1,
    #hjust = 0, 
    #vjust = 0
  #) +
  #geom_segment(
    #data = td_filter(!isTip & node == 293), 
    #aes(x = x + 0.00019, xend = x + 0.00019, y = y, yend = y + 5),  # Increased yend value for longer bar
    #color = "black"
  #) +
  # Label for clade A.2.3
  geom_nodelab(
    data = td_filter(!isTip & node == 602), 
    aes(label = "A.2.3"), 
    geom = "text", 
    nudge_x = 0.00018, 
    nudge_y = 0,
    hjust = 0, 
    vjust = -6
  ) +
  geom_segment(
    data = td_filter(!isTip & node == 602), 
    aes(x = x + 0.000175, xend = x + 0.000175, y = y - 4, yend = y + 50),  # Increased yend value for longer bar
    color = "black"
  ) +
  # Label for clade B.1.2
  geom_nodelab(
    data = td_filter(!isTip & node == 696), 
    aes(label = "B.1.2"), 
    geom = "text", 
    nudge_x = 0.000062,  # Reduced to move label closer to the y-axis
    nudge_y = 0,
    hjust = 0, 
    vjust = -2
  ) +
  geom_segment(
    data = td_filter(!isTip & node == 696), 
    aes(x = x + 0.00006, xend = x + 0.00006, y = y - 5, yend = y + 17),  # Increased yend value for longer bar
    color = "black"
  ) 

print(finalTreePlot$data)
print(finalTreePlot)

# Add necessary columns
treePlotCollapsed$data$xmin <- NA
treePlotCollapsed$data$xmax <- NA
treePlotCollapsed$data$ymin <- NA
treePlotCollapsed$data$ymax <- NA

# Function to get all descendants of a set of nodes
getDescendants <- function(treeDF, nodes) {
  descendants <- c()
  nodes_to_check <- nodes
  
  while (length(nodes_to_check) > 0) {
    current_node <- nodes_to_check[1]
    nodes_to_check <- nodes_to_check[-1]
    child_nodes <- treeDF$node[treeDF$parent == current_node]
    descendants <- c(descendants, child_nodes)
    nodes_to_check <- c(nodes_to_check, child_nodes)
  }
  
  return(descendants)
}

# Define the range of nodes for Clade IIa
nodes_to_include <- 503:516

# Get all descendants of the specified nodes
descendants_of_clade_IIa <- getDescendants(treeDF, nodes_to_include)

# Find the lowest common ancestor (LCA) of the descendants to highlight the entire group
# and not duplicate shading
lca_clade_IIa <- 503  # Replace with the actual LCA node if known

# Highlighting the tree with adjusted Clade IIa coloring
finalTreePlotHighlight <- treePlotCollapsed %<+% cladeDF +
  # Adjust Clade IIa highlight to avoid overlaps by using LCA
  geom_highlight(
    data = treeDF %>% filter(node == lca_clade_IIa), 
    mapping = aes(node = node, fill = "Clade IIa"), 
    type = "roundrect", 
    alpha = 0.8, 
    to.bottom = TRUE
  ) +
  geom_nodelab(data = td_filter(!isTip & node == 502), aes(label = "Clade IIa", color = "black"), geom = "text", nudge_x = 0.000945, nudge_y = -14, hjust = 0, vjust = 0, fontface = "bold") +
  geom_nodelab(data = td_filter(!isTip & node == 547), aes(label = "Clade IIb"), geom = "text", nudge_x = -0.000285, nudge_y = -25, hjust = 0, vjust = 0) +
  geom_nodelab(data = td_filter(!isTip & node == 561), aes(label = "A.2", color = "A.2"), geom = "text", nudge_x = -0.000055, nudge_y = 1, hjust = -0.5, vjust = 0, fontface = "bold") +
  geom_highlight(data = td_filter(!isTip & node == 561), mapping = aes(node = node, fill = "A.2"), type = "roundrect", alpha = 0.2, to.bottom = TRUE) +
  geom_nodelab(data=td_filter(!isTip & node == 647), aes(label = "B.1"), geom="text", nudge_x = -0.000035, nudge_y = 1, hjust= 0, vjust= 0) +
  geom_nodelab(data = td_filter(!isTip & node == 636), aes(label = "A.1"), geom = "text", nudge_x = -0.000035, nudge_y = 1,  hjust = 0, vjust = 0) +
  geom_nodelab(data = td_filter(!isTip & node == 646), aes(label = "A.1.1"), geom = "text", nudge_x = -0.000055, nudge_y = 1,  hjust = 0, vjust = 0) +
  geom_nodelab(data = td_filter(!isTip & node == 528), aes(label = "A"), geom = "text", nudge_x = -0.0000175, nudge_y = 1, hjust = 0, vjust = 0) +
  geom_nodelab(data = td_filter(!isTip & node == 576), aes(label = "A.2.1", color = "A.2.1"), geom = "text", nudge_x = 0.000205, nudge_y = 0, hjust = 0, vjust = -1.35, fontface = "bold") +
  geom_segment(data = td_filter(!isTip & node == 576), aes(x = x + 0.0002, xend = x + 0.0002, y = y - 2, yend = y + 15, color = "A.2.1")) +
  geom_highlight(data = td_filter(!isTip & node == 576), mapping = aes(node = node, fill = "A.2.1"), type = "roundrect", alpha = 1, to.bottom = TRUE) +
  geom_nodelab(data = td_filter(!isTip & node == 602), aes(label = "A.2.3", color = "A.2.3"), geom = "text", nudge_x = 0.00018, nudge_y = 0, hjust = 0, vjust = -6, fontface = "bold") +
  geom_segment(data = td_filter(!isTip & node == 602), aes(x = x + 0.000175, xend = x + 0.000175, y = y - 4, yend = y + 48, color = "A.2.3")) +
  geom_highlight(data = td_filter(!isTip & node == 602), mapping = aes(node = node, fill = "A.2.3"), type = "roundrect", alpha = 1, to.bottom = TRUE) +
  geom_nodelab(data = td_filter(!isTip & node == 696), aes(label = "B.1.2", color = "B.1.2"), geom = "text", nudge_x = 0.000062, nudge_y = 0, hjust = 0, vjust = -2, fontface = "bold") +
  geom_segment(data = td_filter(!isTip & node == 696), aes(x = x + 0.00006, xend = x + 0.00006, y = y - 5, yend = y + 17, color = "B.1.2")) +
  geom_highlight(data = td_filter(!isTip & node == 696), mapping = aes(node = node, fill = "B.1.2"), type = "roundrect", alpha = 0.8, to.bottom = TRUE) +
  scale_fill_manual(values = c("A.2" = "#d4ac0d", "A.2.1" = "#922b21", "A.2.3" = "#4a235a", "Clade IIa" = "#7fb3d5", "B.1.2" = "#229954")) +
  scale_color_manual(values = c("A.2" = "#d4ac0d", "A.2.1" = "#922b21", "A.2.3" = "#4a235a", "Clade IIa" = "#7fb3d5", "B.1.2" = "#229954")) +
  guides(fill = "none", color = "none")

# Print final plot
print(finalTreePlotHighlight)


####################


####################
## This second plot is an example of how to
## zoom in on a certain clade (specifically B.1.2, in
## this case).

# Takes the raw tree plot from earlier and uses the viewClade function 
# to zoom in on the B.1.2 clade. Then the GH samples in that clade are marked with
# red dots, and the points are labled with text. The GH samples are labeled with red
# text and the non-GH samples are labeled with black text
# Define the node for lineage B.1.2
B12Node <- 696

# Recursive function to get descendant tips
get_descendant_tips <- function(tree, node) {
  children <- which(tree$edge[, 1] == node)
  tips <- tree$edge[children, 2]
  
  # If a child is not a tip, recurse
  tips <- unlist(lapply(tips, function(child) {
    if (child <= Ntip(tree)) {  # Check if the node is a tip
      return(child)
    } else {
      return(get_descendant_tips(tree, child))
    }
  }))
  
  return(tips)
}

# Get descendant tips for node 696
descendant_tips <- get_descendant_tips(tree, B12Node)
descendant_labels <- tree$tip.label[descendant_tips]  # Get tip labels

# Subset the tree to include only the descendants of node 696
B12Subtree <- keep.tip(tree, descendant_labels)

# Fortify the subset tree for plotting
B12Data <- ggtree::fortify(B12Subtree)

# Add labels and colors for Ghana and non-Ghana samples
B12Data <- B12Data %>%
  mutate(
    label_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "black"),
    fill_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "white"),
    fontface = ifelse(label %in% GHSamplePoints$label, "bold", "plain") # Bold for Ghana samples
  )

# Plot the zoomed tree with even more right-side space
B12Zoom <- ggtree(B12Subtree) +
  geom_tippoint(
    data = B12Data %>% filter(isTip), 
    aes(x = x, y = y, color = "black", fill = fill_color), # Add fontface for bold text
    size = 2, shape = 21
  ) +
  geom_tiplab(
    data = B12Data %>% filter(isTip), 
    aes(x = x, y = y, label = label, color = label_color, fontface = fontface),
    size = 2.5, 
    nudge_x = 0.00000025  # Fine-tuned label positioning
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +  # Increase right-side space further
  coord_cartesian(clip = "off") +
  ggtitle("Lineage B.1.2") +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Add padding to center the plot
  )

# Display the plot
print(B12Zoom)



####################

####################
## This third plot is zooming in on clade IIa.

getDescendants <- function(treeDF, nodes) {
  descendants <- c()
  nodes_to_check <- nodes
  
  while (length(nodes_to_check) > 0) {
    current_node <- nodes_to_check[1]
    nodes_to_check <- nodes_to_check[-1]
    child_nodes <- treeDF$node[treeDF$parent == current_node]
    descendants <- c(descendants, child_nodes)
    nodes_to_check <- c(nodes_to_check, child_nodes)
  }
  
  return(descendants)
}

# Define the range of nodes
nodes_to_include <- 503:516

# Get all descendants for the range
descendants_of_range <- getDescendants(treeDF, nodes_to_include)

# Filter to find the tip nodes among descendants
descendant_tips_range <- treeDF %>% filter(isTip & node %in% descendants_of_range)

# Subset the tree to include only the desired clade (503-516)
subtree <- keep.tip(tree, tree$tip.label[descendant_tips_range$node])

# Fortify the subtree to ensure proper x and y coordinates
subtree_data <- ggtree::fortify(subtree)

# Add colors for highlighting specific tips
subtree_data <- subtree_data %>%
  mutate(
    tip_color = ifelse(label == "22-010", "#228B22", "black"),  # Highlight sample "22-010" in red
    label_color = ifelse(label == "22-010", "#228B22", "black"), # Label color for "22-010"
    fontface = ifelse(label %in% GHSamplePoints$label, "bold", "plain") # Bold for Ghana samples
  )

# Plot the subtree with customized colors
CIIaZoom <- ggtree(subtree) +
  geom_tippoint(
    data = subtree_data %>% filter(isTip),  # Use only tip data
    aes(
      x = x, y = y, 
      fill = ifelse(label == "22-010", "#228B22", "white"),  # Solid red for 22-010, white for others
      color = "black"  # Keep black border for all points
    ),
    size = 2, shape = 21  # Shape 21 allows for separate fill and border colors
  ) +
  geom_tiplab(
    data = subtree_data %>% filter(isTip),  # Use only tip data
    aes(
      x = x, y = y, 
      color = ifelse(label == "22-010", "#228B22", "black"),  # Red label for 22-010, black for others
      label = label, fontface = fontface
    ),
    size = 2.5,
    nudge_x = 0.000002  # Adjust as needed for label positioning
  ) +
  scale_fill_identity() +  # Use the exact fill colors specified in the dataset
  scale_color_identity() +  # Use the exact label colors specified in the dataset
  coord_cartesian(clip = "off") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 4, 1, 1), "cm")
  ) +
  ggtitle("Clade IIa")

# Display the plot
print(CIIaZoom)



####################
## This fourth plot is zooming in on Lineage A.2.

# Takes the raw tree plot from earlier and uses the viewClade function 
# to zoom in on the Lineage A.2. Then the GH samples in that clade are marked with
# green dots, and the points are labled with bolded text. # Define the node for lineage A.2
A2Node <- 561

# Recursive function to get descendant tips
get_descendant_tips <- function(tree, node) {
  children <- which(tree$edge[, 1] == node)
  tips <- tree$edge[children, 2]
  
  # If a child is not a tip, recurse
  tips <- unlist(lapply(tips, function(child) {
    if (child <= Ntip(tree)) {  # Check if the node is a tip
      return(child)
    } else {
      return(get_descendant_tips(tree, child))
    }
  }))
  
  return(tips)
}

# Get descendant tips for node 561
descendant_tips <- get_descendant_tips(tree, A2Node)
descendant_labels <- tree$tip.label[descendant_tips]  # Get tip labels

# Subset the tree to include only the descendants of node 561
A2Subtree <- keep.tip(tree, descendant_labels)

# Fortify the subset tree for plotting
A2Data <- ggtree::fortify(A2Subtree)

# Add labels and colors for Ghana and non-Ghana samples
A2Data <- A2Data %>%
  mutate(
    label_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "black"),
    fill_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "white"),
    fontface = ifelse(label %in% GHSamplePoints$label, "bold", "plain") # Bold for Ghana samples
  )

# Plot the zoomed tree for Lineage A.2
A2Zoom <- ggtree(A2Subtree) +
  geom_tippoint(
    data = A2Data %>% filter(isTip), 
    aes(x = x, y = y, color = "black", fill = fill_color), 
    size = 2, shape = 21
  ) +
  geom_tiplab(
    data = A2Data %>% filter(isTip), 
    aes(x = x, y = y, label = label, color = label_color, fontface = fontface),
    size = 2.5, 
    nudge_x = 0.000001
  ) +
  
 # Adjust scales for proper spacing
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +  # Add extra space on the right
  coord_cartesian(clip = "off") +
  
  # Add title and adjust theme
  ggtitle("Lineage A.2") +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Add padding to center the plot
  )

# Display the plot
print(A2Zoom)

####################


####################
## This fifth plot is zooming in on Lineage A.2.1 and A.2.3 samples in closer view.

# Function to create a zoomed-in plot for a given node
create_zoom_plot <- function(tree, node, title, GHSamplePoints) {
  # Recursive function to get descendant tips
  get_descendant_tips <- function(tree, node) {
    children <- which(tree$edge[, 1] == node)
    tips <- tree$edge[children, 2]
    
    # If a child is not a tip, recurse
    tips <- unlist(lapply(tips, function(child) {
      if (child <= Ntip(tree)) {  # Check if the node is a tip
        return(child)
      } else {
        return(get_descendant_tips(tree, child))
      }
    }))
    
    return(tips)
  }
  
  # Get descendant tips for the node
  descendant_tips <- get_descendant_tips(tree, node)
  descendant_labels <- tree$tip.label[descendant_tips]  # Get tip labels
  
  # Subset the tree to include only the descendants of the node
  Subtree <- keep.tip(tree, descendant_labels)
  
  # Fortify the subset tree for plotting
  SubtreeData <- ggtree::fortify(Subtree)
  
  # Add labels and colors for Ghana and non-Ghana samples
  SubtreeData <- SubtreeData %>%
    mutate(
      label_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "black"),
      fill_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "white"),
      fontface = ifelse(label %in% GHSamplePoints$label, "bold", "plain") # Bold for Ghana samples
    )
  
  # Plot the zoomed tree
  ZoomPlot <- ggtree(Subtree) +
    geom_tippoint(
      data = SubtreeData %>% filter(isTip), 
      aes(x = x, y = y, color = "black", fill = fill_color), 
      size = 2, shape = 21
    ) +
    geom_tiplab(
      data = SubtreeData %>% filter(isTip), 
      aes(x = x, y = y, label = label, color = label_color, fontface = fontface),
      size = 2.5, 
      nudge_x = 0.000001
    ) +
    
    # Adjust scales for proper spacing
    scale_fill_identity() +
    scale_color_identity() +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +  # Add extra space on the right
    coord_cartesian(clip = "off") +
    
    # Add title and adjust theme
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.margin = unit(c(1, 1, 1, 1), "cm") # Add padding to center the plot
    )
  
  return(ZoomPlot)
}


# Zoom for A.2.1
A21Zoom <- create_zoom_plot(tree, 576, "Lineage A.2.1", GHSamplePoints)
print(A21Zoom)

# Zoom for A.2.3
A23Zoom <- create_zoom_plot(tree, 602, "Lineage A.2.3", GHSamplePoints)
print(A23Zoom)
####################


####################
## This sixth plot is zooming in on just Lineage A.2 samples in closer view.
# Define the nodes
A2Node <- 561
ExcludeNodes <- c(576, 563, 602)  # Nodes to exclude

# Recursive function to get descendant tips
get_descendant_tips <- function(tree, node) {
  children <- which(tree$edge[, 1] == node)
  tips <- tree$edge[children, 2]
  
  # If a child is not a tip, recurse
  tips <- unlist(lapply(tips, function(child) {
    if (child <= Ntip(tree)) {  # Check if the node is a tip
      return(child)
    } else {
      return(get_descendant_tips(tree, child))
    }
  }))
  
  return(tips)
}

# Get all descendants of node 561
descendant_tips <- get_descendant_tips(tree, A2Node)
descendant_labels <- tree$tip.label[descendant_tips]  # Labels for all descendants of 561

# Get all descendants of the exclusion nodes
exclude_tips <- unlist(lapply(ExcludeNodes, function(node) get_descendant_tips(tree, node)))
exclude_labels <- tree$tip.label[exclude_tips]  # Labels for all descendants of exclusion nodes
print("Excluded labels:")
print(exclude_labels)

# Retain descendants of 561, excluding those under the exclusion nodes
filtered_labels <- setdiff(descendant_labels, exclude_labels)
print("Filtered labels:")
print(filtered_labels)

# Subset the tree to include only the filtered labels
A2Subtree <- keep.tip(tree, filtered_labels)

# Ensure the subset tree is not empty
if (is.null(A2Subtree) || Ntip(A2Subtree) == 0) {
  stop("The subset tree is empty. Check the filtering logic.")
}

# Fortify the subset tree for plotting
A2Data <- as_tibble(ggtree::fortify(A2Subtree))

# Add labels and colors for Ghana and non-Ghana samples
A2Data <- A2Data %>%
  mutate(
    label_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "black"),
    fill_color = ifelse(label %in% GHSamplePoints$label, "#228B22", "white"),
    fontface = ifelse(label %in% GHSamplePoints$label, "bold", "plain") # Bold for Ghana samples
  )

# Plot the zoomed tree for Lineage A.2
A2ZoomFiltered <- ggtree(A2Subtree) +
  geom_tippoint(
    data = A2Data %>% filter(isTip), 
    aes(x = x, y = y, color = "black", fill = fill_color), 
    size = 2, shape = 21
  ) +
  geom_tiplab(
    data = A2Data %>% filter(isTip), 
    aes(x = x, y = y, label = label, color = label_color, fontface = fontface),
    size = 2.5, 
    nudge_x = 0.000001
  ) +
  
  # Adjust scales for proper spacing
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +  # Add extra space on the right
  coord_cartesian(clip = "off") +
  
  # Add title and adjust theme
  ggtitle("Lineage A.2 (Excluding Sublineages)") +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    plot.margin = unit(c(1, 1, 1, 1), "cm") # Add padding to center the plot
  )

# Display the plot
print(A2ZoomFiltered)

####################

## UNUSED Code from testing and development of this script

#p <- ggtree(tree) + geom_treescale() +
#  geom_tippoint(data=td_filter(isTip & !(node %in% NESamplePoints$Node)),size=1.25, shape=21, color="black", fill="white") + 
#  geom_tippoint(data=td_filter(isTip & node %in% NESamplePoints$Node),size=2, shape=21, color="black", fill="#f03e47") + coord_cartesian(clip="off") +
#  #geom_tiplab(data=td_filter(isTip & node %in% NESamplePoints$Node), size=2)
#p
#collapse(p, node=nodesToCollapse)

#p2 <- collapseMultipleNodes(p, nodesToCollapse)
#p2 <- collapseMultipleNodes(p, collapseAllB1Sublins)

#p2

#ggsave(filename="test-save.jpg", plot=p2, height=20, width=12, dpi=300)


#highlightData <- NESampleTreeData %>%
#  distinct(lineageParentNode, lineage)


#p2 + geom_hilight(data=highlightData, aes(node=lineageParentNode, fill=lineage), alpha=0.3)

#labData

#p2 + geom_cladelab(date)

# Test tree + zoom ins
#p3 <- ggtree(tree) +
#  geom_hilight(data=highlightData, aes(node=lineageParentNode, fill=lineage), alpha=0.7) + 
#  scale_fill_brewer(palette = "Dark2") + 
#  coord_cartesian(clip="off")

#p3$layers <- rev(p3$layers)

#p3 + geom_tippoint(data=td_filter(isTip & node %in% NESamplePoints$Node),size=2, shape=21, color="black", fill="white") +
#  geom_tiplab(data=td_filter(isTip & node == 1), size=2.5) + geom_treescale(x=0, y=225)


#nonB1Clades <- setdiff(lineageNodeDF$Node, c(1233, 1232))
#treeCollapsed <- collapseMultipleNodes(ggtree(tree), nonB1Clades)
#b1Plot <- viewClade(treeCollapsed, 1233)
#b1Plot + geom_tippoint(data=td_filter(isTip & node %in% NESamplePoints$Node),size=2, shape=21, color="black", fill="white") +
#  coord_cartesian(clip="off")


# Ensure necessary columns are present
#if (!("xmin" %in% colnames(treePlotCollapsed$data))) {
 # treePlotCollapsed$data$xmin <- NA
#}

#if (!("xmax" %in% colnames(treePlotCollapsed$data))) {
 # treePlotCollapsed$data$xmax <- NA
#}

#if (!("ymin" %in% colnames(treePlotCollapsed$data))) {
 # treePlotCollapsed$data$ymin <- NA
#}

#if (!("ymax" %in% colnames(treePlotCollapsed$data))) {
 # treePlotCollapsed$data$ymax <- NA
#}



# Takes the raw tree plot from earlier and uses the viewClade function 
# to zoom in on the Lineage A sample. Then the GH samples in that clade are marked with
# red dots, and the points are labled with text. The GH samples are labeled with red
# text and the non-GH samples are labeled with black text
#AZoom <- viewClade(rawTreePlot, lineageNodeDF$node[lineageNodeDF$lineage == "A"]) + 
#geom_tippoint(data=td_filter(isTip & !(node %in% GHSamplePoints$node)), size=2, shape=21, color="black", fill="white") +
#geom_tiplab(data=td_filter(isTip & node %in% GHSamplePoints$node), color= "#f03e47", size=2.5, nudge_x = 0.000001) +
#geom_tiplab(data=td_filter(isTip & !(node %in% GHSamplePoints$node)), color= "black", size=2.5, nudge_x = 0.000001) +
#ggtitle("Lineage A") + theme(plot.title = element_text(hjust=0.5))

#AZoom 

# Function to get the closest relatives
#getClosestRelatives <- function(target_sample, dist_matrix, n = 10) {
#target_distances <- dist_matrix[target_sample, ]
#closest_nodes <- sort(target_distances, decreasing = FALSE)[1:n]
#closest_nodes_names <- names(closest_nodes)
#return(closest_nodes_names)
#}

# Find the node for the target sample MPXV-22-013
#target_sample <- "MPXV-22-013"
#target_node <- treeDF$node[treeDF$label == target_sample]
#print(target_node)

# Compute the distance matrix
#dist_matrix <- cophenetic(tree)

# Get the closest 10 relatives of the target node
#closest_relatives <- getClosestRelatives(target_sample, dist_matrix, n = 10)
#print(closest_relatives)

# Function to filter the tree based on selected nodes
#filterTree <- function(tree, selected_nodes) {
#tips_to_keep <- intersect(tree$tip.label, selected_nodes)
#filtered_tree <- drop.tip(tree, setdiff(tree$tip.label, tips_to_keep))
#return(filtered_tree)
#}

# Filter the tree to include only the target sample and its closest relatives
#filtered_tree <- filterTree(tree, closest_relatives)

# Plot the filtered tree
#AZoom <- ggtree(filtered_tree) + 
# Points for non-GH samples (white filled)
#geom_tippoint(aes(subset = !(label %in% GHSampleNames)), size = 2, shape = 21, color = "black", fill = "white") +
# Points for GH samples (red filled)
#geom_tippoint(aes(subset = label %in% GHSampleNames), size = 2, shape = 21, color = "black", fill = "#f03e47") +
# Labels for GH samples (red color)
#geom_tiplab(aes(subset = label %in% GHSampleNames, label = label), color = "#f03e47", size = 2.5, nudge_x = 0.00000015) +
# Labels for non-GH samples (black color)
#geom_tiplab(aes(subset = !(label %in% GHSampleNames), label = label), color = "black", size = 2.5, nudge_x = 0.00000015) +
#ggtitle("Lineage A") + theme(plot.title = element_text(hjust = 0.5))

#print(AZoom)

####################

####################

# Finally, creates a plot grid containing the main tree plot and the 
# zoom next to one another.
# Finally, creates a plot grid containing the main tree plot and the 
# zoom next to one another.

# Combine with B12Zoom
#finalPlot_B12 <- plot_grid(finalTreePlotHighlight, B12Zoom, ncol = 2, rel_widths = c(3, 1), labels = c("A", "B"))
#print(finalPlot_B12)

# Combine with CIIaZoom
#finalPlot_CIIa <- plot_grid(finalTreePlotHighlight, CIIaZoom, ncol = 2, rel_widths = c(2, 1), labels = c("A", "C"))
#print(finalPlot_CIIa)

# Combine with A2Zoom
#finalPlot_A2 <- plot_grid(finalTreePlotHighlight, A2Zoom, ncol = 2, rel_widths = c(3, 1), labels = c("A", "D"))
#print(finalPlot_A2)

# Combine with AZoom
#finalPlot_A <- plot_grid(finalTreePlotHighlight, AZoom, ncol = 2, rel_widths = c(2, 1), labels = c("A", "E"))
#print(finalPlot_A)


# Combine the zoomed plots into a single vertical plot
#zoomed_plots_combined <- plot_grid(B12Zoom, CIIaZoom, A2Zoom, ncol = 1, rel_heights = c(1, 1, 1, 1), labels = c("B", "C", "D"))

# Combine the final tree plot with the combined zoomed plots
#final_large_plot <- plot_grid(finalTreePlotHighlight, zoomed_plots_combined, ncol = 2, rel_widths = c(2, 1), labels = c("A"))

# Display the final large plot
#print(final_large_plot)

# Combine B12Zoom and CIIaZoom horizontally
#top_row <- plot_grid(
#B12Zoom + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)),
#CIIaZoom + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)),
#ncol = 2,
#labels = c("B", "C")
#)

# Adjust the title position of A2Zoom
#A2Zoom <- A2Zoom + ggtitle("Lineage A.2") + theme(plot.title = element_text(hjust = 0.5))

# Combine A2Zoom and any other zoomed plots you have vertically
#bottom_row <- plot_grid(
#A2Zoom + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)),
#ncol = 2,
#labels = c("D")
#)

# Combine the top and bottom rows into a single plot vertically
#zoomed_plots_combined <- plot_grid(
#top_row,
#bottom_row,
#ncol = 2,
#rel_heights = c(1, 1)
#)

# Combine the final tree plot with the combined zoomed plots
#final_large_plot <- plot_grid(
#finalTreePlotHighlight + theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10)),
#zoomed_plots_combined,
#ncol = 2,
#rel_heights = c(1, 1),
#labels = c("A")
#)

# Display the final large plot
#print(final_large_plot)
