/*------------------------------------------------------------------------------
* convkml.c : google earth kml converter
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants -----------------------------------------------------------------*/
#define SIZP     0.2            /* mark size of rover positions */
#define SIZR     0.3            /* mark size of reference position */
#define TINT     60.0           /* time label interval (sec) */

static const char *head1="<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
static const char *head2="<kml xmlns=\"http://earth.google.com/kml/2.1\">";
static const char *mark="http://maps.google.com/mapfiles/kml/pal2/icon18.png";

/* output track --------------------------------------------------------------*/
static void outtrack(FILE *f, const solbuf_t *solbuf, const char *color,
                     int outalt, int outtime)
{
    double pos[3];
    int i;

    fprintf(f,"<Placemark>\n");
    fprintf(f,"<name>Rover Track</name>\n");
    fprintf(f,"<Style>\n");
    fprintf(f,"<LineStyle>\n");
    fprintf(f,"<color>%s</color>\n",color);
    fprintf(f,"</LineStyle>\n");
    fprintf(f,"</Style>\n");
    fprintf(f,"<LineString>\n");
    if (outalt) fprintf(f,"<altitudeMode>absolute</altitudeMode>\n");
    fprintf(f,"<coordinates>\n");
    for (i=0;i<solbuf->n;i++) {
        ecef2pos(solbuf->data[i].rr,pos);
        if      (outalt==0) pos[2]=0.0;
        else if (outalt==2) ;
        fprintf(f,"%13.9f,%12.9f,%5.3f\n",pos[1]*R2D,pos[0]*R2D,pos[2]);
    }
    fprintf(f,"</coordinates>\n");
    fprintf(f,"</LineString>\n");
    fprintf(f,"</Placemark>\n");
}
/* output point --------------------------------------------------------------*/
static void outpoint(FILE *fp, gtime_t time, const double *pos,
                     const char *label, int style, int outalt, int outtime, int stat)
{
    double ep[6],alt=0.0;
    char str[256]="";

    fprintf(fp,"<Placemark>\n");
    if (*label) fprintf(fp,"<name>%s</name>\n",label);

    gtime_t utct=gpst2utc(time);

    utct=timeadd(utct,8*3600.0);
    time2epoch(utct,ep);

    fprintf(fp,"<description><![CDATA[<TABLE><TR><TD width=1000>UTC+8: %4d-%02d-%02d  %02d:%02d:%.3lf"
               "<BR>Lat: %12.8f N <BR>Lon: %12.8f E <BR>Hgt: %12.3f m<BR>From: %s</TD></TR></TABLE>]]></description>\n",
               (int)ep[0],(int)ep[1],(int)ep[2],(int)ep[3],(int)ep[4],ep[5],
               pos[0]*R2D,pos[1]*R2D,pos[2],
               stat==SOLQ_FIX?"FIX":stat==SOLQ_FLOAT?"FLOAT":stat==SOLQ_DGPS?"DGNSS":stat==SOLQ_SINGLE?"SINGLE":"NONE");

    fprintf(fp,"<styleUrl>#P%d</styleUrl>\n",style);
    if (outtime) {
        if      (outtime==2) time=gpst2utc(time);
        else if (outtime==3) time=timeadd(gpst2utc(time),9*3600.0);
        time2epoch(time,ep);
        if (!*label&&fmod(ep[5]+0.005,TINT)<0.01) {
            sprintf(str,"%02.0f:%02.0f",ep[3],ep[4]);
            fprintf(fp,"<name>%s</name>\n",str);
        }
        sprintf(str,"%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%05.2fZ",
                ep[0],ep[1],ep[2],ep[3],ep[4],ep[5]);
        fprintf(fp,"<TimeStamp><when>%s</when></TimeStamp>\n",str);
    }
    fprintf(fp,"<Point>\n");
    if (outalt) {
        fprintf(fp,"<extrude>1</extrude>\n");
        fprintf(fp,"<altitudeMode>absolute</altitudeMode>\n");
        alt=pos[2];
    }
    fprintf(fp,"<coordinates>%13.9f,%12.9f,%5.3f</coordinates>\n",pos[1]*R2D,
            pos[0]*R2D,alt);
    fprintf(fp,"</Point>\n");
    fprintf(fp,"</Placemark>\n");
}
/* GNSS solutions convert to KML----------------------------------------------*/
extern int sol2kml(const char *file, const solbuf_t *solbuf, int tcolor,
                   int pcolor, int outalt, int outtime)
{
    FILE *fp;
    double pos[3];
    int i,qcolor[]={0,1,2,5,4,3,0};
    char *color[]={
            "ffffffff","ff008800","ff00aaff","ff0000ff","ff00ffff","ffff00ff"
    };
    if (!(fp=fopen(file,"w"))) {
        fprintf(stderr,"file open error : %s\n",file);
        return 0;
    }
    fprintf(fp,"%s\n%s\n",head1,head2);
    fprintf(fp,"<Document>\n");
    for (i=0;i<6;i++) {
        fprintf(fp,"<Style id=\"P%d\">\n",i);
        fprintf(fp,"  <IconStyle>\n");
        fprintf(fp,"    <color>%s</color>\n",color[i]);
        fprintf(fp,"    <scale>%.1f</scale>\n",i==0?SIZR:SIZP);
        fprintf(fp,"    <Icon><href>%s</href></Icon>\n",mark);
        fprintf(fp,"  </IconStyle>\n");
        fprintf(fp,"</Style>\n");
    }
    if (tcolor>0) {
        outtrack(fp,solbuf,color[tcolor-1],outalt,outtime);
    }
    if (pcolor>0) {
        fprintf(fp,"<Folder>\n");
        fprintf(fp,"  <name>Rover Position</name>\n");
        for (i=0;i<solbuf->n;i++) {
            ecef2pos(solbuf->data[i].rr,pos);

            int clr=0;
            if (solbuf->data[i].stat==SOLQ_FIX ) clr=1;
            if (solbuf->data[i].stat==SOLQ_DGPS) clr=2;
            outpoint(fp,solbuf->data[i].time,pos,"",qcolor[clr],outalt,outtime,solbuf->data[i].stat);
        }
        fprintf(fp,"</Folder>\n");
    }
    fprintf(fp,"</Document>\n");
    fprintf(fp,"</kml>\n");
    fclose(fp);
    return 1;
}

