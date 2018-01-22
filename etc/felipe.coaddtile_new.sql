-- make sure we delete before we created
drop   table felipe.coaddtile_new purge;
create table felipe.coaddtile_new (
TILENAME                VARCHAR2(50),
RA                      NUMBER(15,10),
DEC                     NUMBER(15,10),
RAC1			NUMBER(15,10),
RAC2			NUMBER(15,10),
RAC3			NUMBER(15,10),
RAC4			NUMBER(15,10),
DECC1			NUMBER(15,10),
DECC2			NUMBER(15,10),
DECC3			NUMBER(15,10),
DECC4			NUMBER(15,10),
URAL                    NUMBER(15,10),
URAU                    NUMBER(15,10),
UDECL                   NUMBER(15,10),
UDECU                   NUMBER(15,10),
CROSSRAZERO             char(1) check (CROSSRAZERO in ('N','Y')),
PIXELSCALE              NUMBER(10,5),
NAXIS1                  NUMBER(6),
NAXIS2                  NUMBER(6),
CRPIX1                  NUMBER(6),
CRPIX2                  NUMBER(6),
CRVAL1                  NUMBER(15,10),
CRVAL2                  NUMBER(15,10),
CD1_1                   NUMBER(15,10),
CD1_2                   NUMBER(15,10),
CD2_1                   NUMBER(15,10),
CD2_2                   NUMBER(15,10),
constraint coaddtile_new PRIMARY KEY (tilename)
);
-- Add description of columns
comment on column felipe.coaddtile_new.DEC    is 'RA  center of DES tile (deg)'; 
comment on column felipe.coaddtile_new.RA     is 'DEC center of DES file (deg)';
comment on column felipe.coaddtile_new.RAC1   is 'Corner 1 RA of DES tile (deg)';
comment on column felipe.coaddtile_new.RAC2   is 'Corner 2 RA of DES tile (deg)';
comment on column felipe.coaddtile_new.RAC3   is 'Corner 3 RA of DES tile (deg)';
comment on column felipe.coaddtile_new.RAC4   is 'Corner 4 RA of DES tile (deg)';
comment on column felipe.coaddtile_new.DECC1  is 'Corner 1 DEC of DES tile (deg)';
comment on column felipe.coaddtile_new.DECC2  is 'Corner 2 DEC of DES tile (deg)';
comment on column felipe.coaddtile_new.DECC3  is 'Corner 3 DEC of DES tile (deg)';
comment on column felipe.coaddtile_new.DECC4  is 'Corner 4 DEC of DES tile (deg)';
comment on column felipe.coaddtile_new.CROSSRAZERO  is 'DES tile crosses RA=0 [Y/N]';
