import pypistats as pis
from datetime import date, timedelta
import sqlite3

if __name__ == '__main__':
    dbfile = 'prody_stats.db'
    package = 'prody'
    mirror = False
    base_date = "2011-10-01"
    base_number = 2135029

    # set up time elapse from yesterday to today
    today = date.today()
    str_today = today.strftime('%Y-%m-%d')

    yesterday = date.today() - timedelta(days=1)
    str_yesterday = yesterday.strftime('%Y-%m-%d')

    # retrieve download counts for yesterday
    ret = pis.overall(package, mirrors=mirror, 
                               start_date=str_yesterday, 
                               end_date=str_today, 
                               format='numpy')
    n_downloads = ret[0][-1]
    
    # connecting to the database  
    connection = sqlite3.connect(dbfile) 
    crsr = connection.cursor() 
    
    # create the main table (if not already exists) in the database 
    sql_cmd = 'CREATE TABLE IF NOT EXISTS main (date DATE PRIMARY KEY, downloads INTEGER);'
    crsr.execute(sql_cmd)
    
    # insert the base number 
    sql_cmd = 'INSERT OR REPLACE INTO main VALUES ("{0}", {1});'.format(base_date, base_number)
    crsr.execute(sql_cmd)

    # insert or update the data in the table 
    sql_cmd = 'INSERT OR REPLACE INTO main VALUES ("{0}", {1});'.format(str_yesterday, n_downloads)
    crsr.execute(sql_cmd)

    # summation
    sql_cmd = 'SELECT SUM(downloads) FROM main;'
    crsr.execute(sql_cmd)
    ret = crsr.fetchone()
    print(ret[0])
    
    # save and close the database
    connection.commit() 
    connection.close() 

