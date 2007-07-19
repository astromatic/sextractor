<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
        <!ENTITY nbsp "&#160;">
        <!ENTITY deg "&#176;">
        <!ENTITY amin "&#180;">
        <!ENTITY asec "&#168;">
        ]>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- *********************** Global XSL template ************************** -->
 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <html>

<!-- HTML head  -->

   <head>

<!-- javascript -->

<!--  <script type="text/javascript" language="javascript"> -->

    <script src="http://terapix.iap.fr/cplt/xsl/sorttable.js"/>

    <style type="text/css">
     p.sansserif {font-family: sans-serif}
     body {background-color: white}
     mono {font-family: monospace}
     elen {font-family: monospace; font-size: 100%; font-weight: bold; color: green }
     elep {font-family: monospace; font-size: 100%; font-weight: bold; color: red }
     el {font-family: monospace; font-size: 100%; color: black}
     a {text-decoration: none}
     table.sortable a.sortheader
      {
      background-color:#FFEECC;
      color: black;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }
     table.sortable span.sortarrow
      {
      color: black;
      font-weight: bold;
      text-decoration: none;
      }
     table.sortable a.sortheader.sub
      {
      vertical-align: sub;
      }
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </head>

<!-- HTML body -->

    <BODY>
     <TABLE BORDER="0" CELLPADDING="0" CELLSPACING="0" WIDTH="100%">
      <TR>
       <TD ALIGN="LEFT">
        <TABLE BORDER="0">
         <TR>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix"/>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo"/>
          </TD>
          <TD ALIGN="CENTER">
           <FONT color="#669933">
            <B> Processing summary</B>
           </FONT>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixPicture.gif" ALT="Terapix banner"/>
          </TD>
         </TR>
        </TABLE>
       </TD>
      </TR>
      <TR>
       <TD>
        <TABLE BORDER="0" WIDTH="100%" BGCOLOR="#000000">
         <TR>
          <TH BGCOLOR="#000000" ALIGN="LEFT"><FONT SIZE="-1" COLOR="#FFFFFF"> Home > Tools > Data reduction</FONT></TH>
         </TR>
        </TABLE>
       </TD>
      </TR>
     </TABLE>
    <xsl:call-template name="VOTable"/>
   </BODY>
  </html>
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
  <p>
  <sans-serif><i>click to expand or hide tables</i></sans-serif>
  </p>
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
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('extdata')">
     Summary Table on Output Catalog
    </BUTTON>
    <TABLE class="sortable" id="extdata" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Extension</TH>
      <TH BGCOLOR="#FFEECC">Date</TH>
      <TH BGCOLOR="#FFEECC">Time</TH>
      <TH BGCOLOR="#FFEECC">Duration</TH>
      <TH BGCOLOR="#FFEECC">Detected Source Number</TH>
      <TH BGCOLOR="#FFEECC">Sextracted Source Number</TH>
      <TH BGCOLOR="#FFEECC">Image ID</TH>
      <TH BGCOLOR="#FFEECC">Mean Background</TH>
      <TH BGCOLOR="#FFEECC">Sigma Background</TH>
      <TH BGCOLOR="#FFEECC">Detection Threshold</TH>
      <TH BGCOLOR="#FFEECC">Weight Scaling</TH>
      <TH BGCOLOR="#FFEECC">Pixel Scale</TH>
      <TH BGCOLOR="#FFEECC">Epoch</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$extension]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$duration]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ndet]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nsex]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$imid]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$backmean],'#####0.00')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$backsig],'###0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$thresh],'#####0.00')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$wscale],'#####0.00')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$pixscale],'#0.000')"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
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
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('config')">
     Configuration File: <xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/>
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Config Parameter</TH>
      <TH BGCOLOR="#FFEECC">Value</TH>
     </TR>
     <tr BGCOLOR="#EEEEEE">
      <td><el>CATALOG_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Catalog_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>CATALOG_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Catalog_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PARAMETERS_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Parameters_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>DETECT_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Detect_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>DETECT_MINAREA</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Detect_MinArea']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>THRESH_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Thresh_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>DETECT_THRESH</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Detect_Thresh']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ANALYSIS_THRESH</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Analysis_Thresh']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>FILTER</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Filter']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>FILTER_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Filter_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>DEBLEND_NTHRESH</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Deblend_NThresh']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>DEBLEND_MINCONT</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Deblend_MinCont']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>CLEAN</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Clean']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>CLEAN_PARAM</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Clean_Param']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>MASK_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Mask_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>WEIGHT_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Weight_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>WEIGHT_THRESH</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Weight_Thresh']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>WEIGHT_IMAGE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Weight_Image']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>WEIGHT_GAIN</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Weight_Gain']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>FLAG_IMAGE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Flag_Image']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>FLAG_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Flag_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PHOT_APERTURES</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Phot_Apertures']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PHOT_AUTOPARAMS</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Phot_AutoParams']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PHOT_PETROPARAMS</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Phot_PetroParams']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PHOT_AUTOAPERS</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Phot_AutoApers']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PHOT_FLUXFRAC</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Phot_FluxFrac']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>SATUR_LEVEL</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Satur_Level']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>MAG_ZEROPOINT</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Mag_ZeroPoint']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>MAG_GAMMA</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Mag_Gamma']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>GAIN</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Gain']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PIXEL_SCALE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Pixel_Scale']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>SEEING_FWHM</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Seeing_FWHM']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>STARNNW_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='StarNNW_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>BACK_SIZE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Back_Size']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>BACK_FILTERSIZE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Back_FilterSize']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>BACKPHOTO_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='BackPhoto_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>BACKPHOTO_THICK</el></td>
      <td><el><xsl:value-of select="PARAM[@name='BackPhoto_Thick']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>CHECKIMAGE_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='CheckImage_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>CHECKIMAGE_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='CheckImage_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>MEMORY_OBJSTACK</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Memory_ObjStack']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>MEMORY_PIXSTACK</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Memory_PixStack']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>MEMORY_BUFSIZE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Memory_BufSize']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ASSOC_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Assoc_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ASSOC_DATA</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Assoc_Data']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ASSOC_PARAMS</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Assoc_Params']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ASSOC_RADIUS</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Assoc_Radius']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ASSOC_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Assoc_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>ASSOCSELEC_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='AssocSelec_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>VERBOSE_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='Verbose_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>FITS_UNSIGNED</el></td>
      <td><el><xsl:value-of select="PARAM[@name='FITS_Unsigned']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PSF_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='PSF_Name']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PSF_NMAX</el></td>
      <td><el><xsl:value-of select="PARAM[@name='PSF_NMax']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>PSFDISPLAY_TYPE</el></td>
      <td><el><xsl:value-of select="PARAM[@name='PSFDisplay_Type']/@value"/></el></td>
     </tr>
     <tr BGCOLOR="#EEEEEE">
      <td><el>SOM_NAME</el></td>
      <td><el><xsl:value-of select="PARAM[@name='SOM_Name']/@value"/></el></td>
     </tr>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('commandline')">
     Command Line
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD BGCOLOR="#FFEECC" style="font-size: 80%;"><el>Command Line: <xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el></TD>
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
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('warnings')">
     Warnings (limited to the last 100)
    </BUTTON>
    <TABLE id="warnings" style="display: none">
     <TR style="font-size: 80%;">
      <TH BGCOLOR="#FFEECC">Date</TH>
      <TH BGCOLOR="#FFEECC">Time</TH>
      <TH BGCOLOR="#FFEECC">Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
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
