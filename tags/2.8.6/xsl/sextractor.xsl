<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
        <!ENTITY nbsp "&#160;">
        <!ENTITY deg "&#176;">
        <!ENTITY amin "&#180;">
        <!ENTITY asec "&#168;">
        <!ENTITY darr "&#8595;">
        ]>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<!-- ****************** Global XSL template for SExtractor ****************
     (C) E.Bertin and C.Marmo IAP/CNRS/UPMC 2005-2009
     ********************************************************************** -->

 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <HTML>
   <HEAD>
    <link rel="shortcut icon" type="image/x-icon" href="http://astromatic.iap.fr/xsl/favicon.ico" />
    <script type="text/javascript" src="http://astromatic.iap.fr/xsl/sorttable.js"/>

    <style type="text/css">
     p {
      font-family: sans-serif;
      }
     p.italic {font-style: italic}
     body {
      margin: 10px;
      background-color: #e0e0e0;
      background-image: url("http://astromatic.iap.fr/xsl/body_bg.jpg");
      background-repeat: repeat-x;
      background-position: top;
      min-width:662px;
      }
     mono {font-family: monospace}
     elmin {color: green }
     elmax {color: red }
     el {
      font-family: monospace;
      font-size: 100%;
      color: black;
      }
     elm {
      font-family: monospace;
      font-size: 67%;
      white-space: nowrap;
      }
     elh {font-family: monospace;
      font-size: 120%;
      }
     elhi {
      font-family: monospace;
      font-size: 100%;
      white-space: nowrap;
      font-weight: normal;
      font-style: italic;
      }
     a {text-decoration: none; font-style: bold; color: #476674}
     a:hover {text-decoration: underline;}
     #header {
      padding: 5px;
      min-width: 662px;
      background-image: url("http://astromatic.iap.fr/xsl/astromaticleft.png");
      background-repeat: repeat-x;
      background-position: left top;
      text-align: left;
      font-size: 1.2em;
      margin: 0 0 30px 0;
      color:#d3e7f0;
      font-weight: bold;
      }
     th {
      background-color:#d3e7f0;
      border-top: 1px solid white;
      border-left: 1px solid white;
      border-right: 1px solid #476674;
      border-bottom: 1px solid #476674;
      -moz-border-radius: 3px;
      -khtml-border-radius: 3px;
      -webkit-border-radius: 3px;
      border-radius: 3px;
      padding: 2px;
      line-height: 12px;
      }
     td {
      background-color:#f2f4f4;
      padding-left: 2px;
      padding-right: 2px;
      }
     table.sortable {
      border-top: 1px solid #476674;
      border-left: 1px solid #476674;
      border-right: 1px solid white;
      border-bottom: 1px solid white;
      -moz-border-radius: 3px;
      -khtml-border-radius: 3px;
      -webkit-border-radius: 3px;
      border-radius: 3px;
      }
     table.sortable a.sortheader {
      background-color:#d3e7f0;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }

     table.sortable span.sortarrow {
      color: black;
      font-weight: bold;
      text-decoration: blink;
      }
     table.sortable a.sortheader.sub {vertical-align: sub}
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </HEAD>
    <BODY>
     <div id="header">
      <a href="/"><img style="vertical-align: middle; border:0px" src="http://astromatic.iap.fr/xsl/astromatic.png" title="Astromatic home" alt="Astromatic.net" /></a>  Processing summary
     </div>
     <xsl:call-template name="VOTable"/>
   </BODY>
  </HTML>
 </xsl:template>

<!-- **************** Generic XSL template for VOTables ****************** -->
 <xsl:template name="VOTable">
  <xsl:for-each select="/VOTABLE">
   <xsl:call-template name="Resource"/>
  </xsl:for-each>
 </xsl:template>

<!-- *************** Generic XSL template for Resources ****************** -->
 <xsl:template name="Resource">
  <xsl:for-each select="RESOURCE">
   <xsl:choose>
    <xsl:when test="@ID='SExtractor'">
     <xsl:call-template name="sextractor"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>

<!-- ******************* XSL template for SExtractor ********************* -->
 <xsl:template name="sextractor">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="/VOTABLE/RESOURCE[@ID='SExtractor']/TABLE[@ID='Source_List']">
    <xsl:call-template name="sources"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Extension_Data']">
    <xsl:call-template name="extdata"/>
   </xsl:for-each>
   <xsl:for-each select="RESOURCE[@ID='Config']">
    <xsl:call-template name="Config"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Warnings']">
    <xsl:call-template name="Warnings"/>
   </xsl:for-each>
  </xsl:for-each>
 </xsl:template>

<!-- ************* Generic XSL RunInfo template for MetaData ************* -->
 <xsl:template name="RunInfo">
  <p>
<!-- Software name, version, date, time and number of threads -->
   <a>
    <xsl:attribute name="href">
     <xsl:value-of select="PARAM[@name='Soft_URL']/@value"/>
    </xsl:attribute>
    <b>
     <xsl:value-of select="PARAM[@name='Software']/@value"/>&nbsp;<xsl:value-of select="PARAM[@name='Version']/@value"/>
    </b>
   </a>
   started on
   <b><xsl:value-of select="PARAM[@name='Date']/@value"/></b>
   at
   <b><xsl:value-of select="PARAM[@name='Time']/@value"/></b>
   with
   <b><xsl:value-of select="PARAM[@name='NThreads']/@value"/></b>
   thread<xsl:if test="PARAM[@name='NThreads']/@value &gt; 1">s</xsl:if>

<!-- Run time -->
   <xsl:variable name="duration" select="PARAM[@name='Duration']/@value"/>
   (run time:
    <b>
     <xsl:choose> 
      <xsl:when test="$duration &gt; 3600.0">
       <xsl:value-of
        select='concat(string(floor($duration div 3600)),
        " h ", format-number(floor(($duration div 60) mod 60.0), "00"),
        " min")'/>
      </xsl:when>
      <xsl:otherwise>
       <xsl:choose>
        <xsl:when test="$duration &gt; 60.0">
         <xsl:value-of
          select='concat(format-number(floor($duration div 60),"##"),
          " min ", format-number(floor($duration mod 60.0), "00")," s")'/>
        </xsl:when>
        <xsl:otherwise>
         <xsl:value-of select='concat(string($duration), " s")'/>
        </xsl:otherwise>
       </xsl:choose>
      </xsl:otherwise>
     </xsl:choose>
    </b>)
    <br />
   by user <b><xsl:value-of select="PARAM[@name='User']/@value"/></b>
   from <b><xsl:value-of select="PARAM[@name='Host']/@value"/></b>
   in <b><mono><xsl:value-of select="PARAM[@name='Path']/@value"/></mono></b>
  </p>
  <p>
   <b style="color: red"><xsl:if test="PARAM[@name='Error_Msg']/@value &gt; 0">
    An Error occured!!! </xsl:if>
   <xsl:value-of select="PARAM[@name='Error_Msg']/@value"/></b>
  </p>
 </xsl:template>

<!-- ******************** XSL template for Source List ******************** -->
 <xsl:template name="sources">
  <xsl:choose> 
   <xsl:when test="DATA/TABLEDATA">
    <p>
     <BUTTON type="button" onclick="showhideTable('sources')" title="click to expand">
     Source List&nbsp;&darr;
     </BUTTON>
     <TABLE id="sources" class="sortable" style="display: none">
      <TR>
       <xsl:for-each select="FIELD">
        <TH align="center"><xsl:attribute name="title"><xsl:value-of select="DESCRIPTION"/></xsl:attribute>
         <elh><xsl:value-of select="@name"/></elh>
         <BR />
         <elhi>
          <xsl:value-of select="@unit"/>
          <xsl:if test="@unit = ''">-</xsl:if>
         </elhi>
        </TH>
       </xsl:for-each>
      </TR>
      <xsl:for-each select="DATA/TABLEDATA">
       <xsl:for-each select="TR">
        <tr>
         <xsl:for-each select="TD">
          <td align="right" >
           <el><xsl:value-of select="self::TD"/></el>
          </td>
         </xsl:for-each>
        </tr>
       </xsl:for-each>
      </xsl:for-each>
     </TABLE>
    </p>
   </xsl:when>
   <xsl:otherwise>
    <p>
     <BUTTON type="button" onclick="showhideTable('catparam')" title="click to expand">
     Parameter List&nbsp;&darr;
     </BUTTON>
     <TABLE id="catparam" class="sortable" style="display: none">
      <TR>
       <TH  align="center">Parameter Name</TH>
       <TH  align="center">Description</TH>
      </TR>
      <xsl:for-each select="FIELD">
       <tr>
        <td align="left" >
         <el><xsl:value-of select="@name"/></el>
        </td>
        <td align="left" >
         <el><xsl:value-of select="DESCRIPTION"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </TABLE>
    </p>
   </xsl:otherwise> 
  </xsl:choose>  
 </xsl:template>

<!-- ******************* XSL template for Extension Data ****************** -->
  <xsl:template name="extdata">
   <xsl:variable name="extension" select="count(FIELD[@name='Extension']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="duration" select="count(FIELD[@name='Duration']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ndet" select="count(FIELD[@name='NDetect']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nsex" select="count(FIELD[@name='NSextracted']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="imid" select="count(FIELD[@name='Image_Ident']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backmean" select="count(FIELD[@name='Background_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="backsig" select="count(FIELD[@name='Background_StDev']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="thresh" select="count(FIELD[@name='Threshold']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="wscale" select="count(FIELD[@name='Weight_Scaling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="pixscale" select="count(FIELD[@name='Pixel_Scale']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="epoch" select="count(FIELD[@name='Epoch']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" onclick="showhideTable('extdata')" title="click to expand">
     Summary Table on Output Catalog&nbsp;&darr;
    </BUTTON>
    <TABLE class="sortable" id="extdata" style="display: none">
     <TR>
      <TH >Extension</TH>
      <TH >Date</TH>
      <TH >Time</TH>
      <TH >Duration</TH>
      <TH >Detected Source Number</TH>
      <TH >Sextracted Source Number</TH>
      <TH >Image ID</TH>
      <TH >Mean Background</TH>
      <TH >Sigma Background</TH>
      <TH >Detection Threshold</TH>
      <TH >Weight Scaling</TH>
      <TH >Pixel Scale</TH>
      <TH >Epoch</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" >
         <el><xsl:value-of select="TD[$extension]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$duration]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$ndet]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$nsex]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$imid]"/></el>
        </td>
        <td align="right" >
         <el><xsl:value-of select="format-number(TD[$backmean],'#####0.00')"/></el>
        </td>
        <td align="right" >
         <el><xsl:value-of select="format-number(TD[$backsig],'###0.000')"/></el>
        </td>
        <td align="right" >
         <el><xsl:value-of select="format-number(TD[$thresh],'#####0.00')"/></el>
        </td>
        <td align="right" >
         <el><xsl:value-of select="format-number(TD[$wscale],'#####0.00')"/></el>
        </td>
        <td align="right" >
         <el><xsl:value-of select="format-number(TD[$pixscale],'#0.000')"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$epoch]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

<!-- ******************** XSL template for Config File ******************** -->
  <xsl:template name="Config">
   <p>
    <BUTTON type="button" onclick="showhideTable('config')" title="click to expand">
     Configuration File:
     <B><xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/></B>
     &darr;
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH >Config Parameter</TH>
      <TH >Value</TH>
     </TR>
     <xsl:for-each select="PARAM[position()>2]">
      <tr >
       <td><el><xsl:value-of select="@name"/></el></td>
       <td><el><xsl:value-of select="@value"/></el></td>
      </tr>
     </xsl:for-each>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" onclick="showhideTable('commandline')" title="click to expand">
     Command Line&nbsp;&darr;
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD  style="font-size: 80%;">
       <el><xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el>
      </TD>
     </TR>
    </TABLE>
   </p>
  </xsl:template>

<!-- ********************** XSL template for Warnings ********************** -->
  <xsl:template name="Warnings">
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="msg" select="count(FIELD[@name='Msg']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" onclick="showhideTable('warnings')" title="click to expand">
     Warnings (limited to the last 100)&nbsp;&darr;
    </BUTTON>
    <TABLE id="warnings" class="sortable" style="display: none">
     <TR style="font-size: 80%;">
      <TH >Date</TH>
      <TH >Time</TH>
      <TH >Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  >
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td >
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" >
         <el><xsl:value-of select="TD[$msg]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>

 <xsl:template name="Rest">
</xsl:template>
</xsl:stylesheet>
