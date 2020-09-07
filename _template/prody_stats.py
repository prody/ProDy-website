import pypistats as pis
from datetime import date, timedelta
from mysql import connector as sql

if __name__ == '__main__':
    host = 'localhost'
    user = 'prody'
    passwd = "I'm Protein Dynamics"
    package = 'prody'
    database = 'prody'
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
    connection = sql.connect(host=host, user=user, passwd=passwd, database=database) 
    crsr = connection.cursor() 
    
    # create the main table (if not already exists) in the database 
    sql_cmd = 'CREATE TABLE IF NOT EXISTS downloads (date DATE PRIMARY KEY, number INTEGER);'
    crsr.execute(sql_cmd)
    
    # insert the base number 
    sql_cmd = 'REPLACE INTO downloads VALUES ("{0}", {1});'.format(base_date, base_number)
    crsr.execute(sql_cmd)

    # insert or update the data in the table 
    sql_cmd = 'REPLACE INTO downloads VALUES ("{0}", {1});'.format(str_yesterday, n_downloads)
    crsr.execute(sql_cmd)

    # summation
    sql_cmd = 'SELECT SUM(number) FROM downloads;'
    crsr.execute(sql_cmd)
    ret = crsr.fetchone()
    print(ret[0])
    
    # save and close the database
    connection.commit() 
    connection.close() 

