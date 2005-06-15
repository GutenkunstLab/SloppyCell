<?xml version="1.0"?>
<!-- 

FILE : sbml_l2v1_todot.xsl

CREATED : January 2004

LAST MODIFIED : 26th April 2005 by Joanne Matthews (j.4.matthews@herts.ac.uk)

AUTHOR : Joanne Matthews(j.4.matthews@herts.ac.uk)
         Biocomputation Research Group
         Science and Technology Research Centre
         University of Hertfordshire, UK

DESCRIPTION : This stylesheet converts an sbml level 2 version 1model into a dot file
CHANGES :
26/04/2005 - JM - Enclosed species names in quotes because label truncates in dotty if hyphen exists in name..

Notes
	
 -->

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.1" xmlns:sbml="http://www.sbml.org/sbml/level2" exclude-result-prefixes="sbml">
<xsl:output method="text"/>
<xsl:strip-space elements="*"/>
<xsl:template match="/">
<xsl:apply-templates/>}
</xsl:template>
<xsl:template match="sbml:model">
digraph "<xsl:choose>
	<xsl:when test="@name"><xsl:value-of select="@name"/>"{
	</xsl:when>
	<xsl:otherwise><xsl:value-of select="@id"/>"{
	</xsl:otherwise>
	</xsl:choose>
<xsl:apply-templates select="sbml:listOfSpecies"/>
<xsl:apply-templates select="sbml:listOfReactions"/>

</xsl:template>

<xsl:template match="sbml:species">
<!--<xsl:template match="sbml:species[@boundaryCondition='true']">-->
<xsl:choose>
	<xsl:when test='@boundaryCondition="true"'>
		<xsl:choose>
			<xsl:when test="@name">
				"<xsl:value-of select="@name"/>"[color=blue]
			</xsl:when>
			<xsl:otherwise>
				"<xsl:value-of select="@id"/>"[color=blue]
			</xsl:otherwise>
		</xsl:choose>
	</xsl:when>
	<xsl:otherwise>
		<xsl:choose>
			<xsl:when test="@name">
				"<xsl:value-of select="@name"/>"[color=black]
			</xsl:when>
			<xsl:otherwise>
				"<xsl:value-of select="@id"/>"[color=black]
			</xsl:otherwise>
		</xsl:choose>
	</xsl:otherwise>
</xsl:choose>
</xsl:template>


<xsl:template match="sbml:reaction">
<xsl:param name="theReaction">
<xsl:choose>
	<xsl:when test="@name">
		<xsl:if test='@name=""'>
			<xsl:value-of select="@id"/>
		</xsl:if>
<xsl:value-of select="@name"/>
	</xsl:when>
<xsl:otherwise>
<xsl:value-of select="@id"/>
</xsl:otherwise>
</xsl:choose>
</xsl:param>
<!--create a rectangular node for the reaction
By default a reaction is reversible. Reversible reactions are red-->
"<xsl:value-of select="$theReaction"/>"[shape=rectangle]
<xsl:choose>
	<xsl:when test="@reversible='false'">
	[color=black]
	</xsl:when>
	<xsl:otherwise>
		[color=red]
	</xsl:otherwise>
</xsl:choose>
<!--
<xsl:if test="@reversible='true'">
[color=red]
</xsl:if>
-->

<xsl:apply-templates select="sbml:listOfReactants"/>
<xsl:apply-templates select="sbml:listOfProducts"/>
<xsl:apply-templates select="sbml:listOfModifiers"/>
</xsl:template>

<xsl:template match="sbml:listOfReactants/sbml:speciesReference">
<!--Get the reactant -->
<xsl:variable name="specName" select="ancestor::sbml:model/sbml:listOfSpecies/sbml:species[@id=current()/@species]"/>

<xsl:choose>
	<xsl:when test="$specName/@name">"<xsl:value-of select="$specName/@name"/>"
	</xsl:when>
<xsl:otherwise>
"<xsl:value-of select="$specName/@id"/>"
</xsl:otherwise>
</xsl:choose>
->
<xsl:choose>
	<xsl:when test="ancestor::sbml:reaction/@name">"<xsl:if test='ancestor::sbml:reaction/@name=""'>
			<xsl:value-of select="ancestor::sbml:reaction/@id"/>
		</xsl:if>
<xsl:value-of select="ancestor::sbml:reaction/@name"/>";
	</xsl:when>
	<xsl:otherwise>
	<xsl:value-of select="ancestor::sbml:reaction/@id"/>;
	</xsl:otherwise>
</xsl:choose>
</xsl:template>

<xsl:template match="sbml:listOfProducts/sbml:speciesReference">
<xsl:choose>
	<xsl:when test="ancestor::sbml:reaction/@name">"<xsl:if test='ancestor::sbml:reaction/@name=""'>
			<xsl:value-of select="ancestor::sbml:reaction/@id"/>
		</xsl:if>
<xsl:value-of select="ancestor::sbml:reaction/@name"/>"->
	</xsl:when>
	<xsl:otherwise>
	<xsl:value-of select="ancestor::sbml:reaction/@id"/>->
	</xsl:otherwise>
</xsl:choose>
<!-- get the product -->
<xsl:variable name="specName" select="ancestor::sbml:model/sbml:listOfSpecies/sbml:species[@id=current()/@species]"/>

<xsl:choose>
	<xsl:when test="$specName/@name">"<xsl:value-of select="$specName/@name"/>";
	</xsl:when>
<xsl:otherwise>
"<xsl:value-of select="$specName/@id"/>";
</xsl:otherwise>
</xsl:choose>


</xsl:template>


<xsl:template match="sbml:modifierSpeciesReference">
<!--<xsl:value-of select="@species"/>->-->

<xsl:variable name="specName" select="ancestor::sbml:model/sbml:listOfSpecies/sbml:species[@id=current()/@species]"/>

<xsl:choose>
	<xsl:when test="$specName/@name">"<xsl:value-of select="$specName/@name"/>"
	</xsl:when>
<xsl:otherwise>
"<xsl:value-of select="$specName/@id"/>"
</xsl:otherwise>
</xsl:choose>
->

<xsl:choose>
	<xsl:when test="ancestor::sbml:reaction/@name">"<xsl:if test='ancestor::sbml:reaction/@name=""'>
			<xsl:value-of select="ancestor::sbml:reaction/@id"/>
		</xsl:if>
<xsl:value-of select="ancestor::sbml:reaction/@name"/>"
	</xsl:when>
	<xsl:otherwise>
	"<xsl:value-of select="ancestor::sbml:reaction/@id"/>"
	</xsl:otherwise>
</xsl:choose>
[style=dotted];
</xsl:template>


</xsl:stylesheet>



