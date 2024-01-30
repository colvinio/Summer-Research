import datetime
import mysql.connector as sql

def main():
    print(datetime.datetime.now().time())

    db = sql.connect(host="ip-173-16.main.oberlin.edu", user="colvin", password="c0LV1n", database="sdss")
    c = db.cursor()

    #Get all objids where there are at least 3 metallicites, currently the limit has been removed
    c.execute(f'''select distinct(r_kpc) from dr17_spaxels_uber limit 10;''')
            
    galaxiesIDS = c.fetchall()
    db.close()

    print(galaxiesIDS)

    exit()

if __name__ == "__main__":
    main()