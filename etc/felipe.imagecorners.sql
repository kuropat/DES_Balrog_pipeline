-- make sure we delete before we created
drop   table felipe.imagecorners purge;
create table felipe.imagecorners (
ID                      NUMBER(22,10),
IMAGENAME               VARCHAR2(100),
RA  			NUMBER(22,8),
DEC  			NUMBER(22,8),
RAC1 			NUMBER(22,8),
RAC2 			NUMBER(22,8),
RAC3 			NUMBER(22,8),
RAC4 			NUMBER(22,8),
DECC1 			NUMBER(22,8),
DECC2 			NUMBER(22,8),
DECC3 			NUMBER(22,8),
DECC4 			NUMBER(22,8),
constraint imagecorners PRIMARY KEY (ID)
);
-- Add description of columns here
