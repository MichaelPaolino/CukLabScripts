% ODBC data source name for raman.accdb
dbDSN = 'Raman Access';

% Establish ODBC connection to data source
conn = database(dbDSN, '','');

% assign foreign key values
fkCols = genFKStruct('UserID','Users','ID','FirstName','Ilya');

procdTableGS = myFSRS(inds,3).dbStageSave(conn,'FSRSProcd','prefix','STO GS','fkStruct',fkCols,'update',356:359);
procdTableES = myFSRS(inds,2).dbStageSave(conn,'FSRSProcd','prefix','STO ES','fkStruct',fkCols,'update',360:363);
procdTableESmGS = myFSRS(inds,1).dbStageSave(conn,'FSRSProcd','prefix','STO ES-GS','fkStruct',fkCols,'update',364:367);
procdTableTR = myTR(inds,1).dbStageSave(conn,'FSRSProcd','prefix','STO TR','fkStruct',fkCols,'update',[368:371]);

% Close the ODBC connection
% conn.close