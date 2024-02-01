import mysql.connector as sql
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pickle

#galaxy LIMIT has been REMOVED

def main():
    print(datetime.datetime.now().time())
    print()

    metal_types = ['D16_smc', 'PP04_O3N2_smc', 'PP04_N2_smc', 'M13_N2_smc', 'M13_O3N2_smc', 'C17_O3N2_smc', 'C17_N2_smc', 'D16_mw', 'PP04_O3N2_mw', 'PP04_N2_mw', 'M13_N2_mw', 'M13_O3N2_mw', 'C17_O3N2_mw', 'C17_N2_mw'] # Z94_smc, 'M91_smc', 'KK04_smc', 'KE08_smc', 'Z94_mw', 'M91_mw', 'KK04_mw', 'KE08_mw', 

    for type in metal_types:
        
        #Resetting universe data dictionary for each metallicity
        universe_data = {}
        
        print(f"Pulling objids for {type}...")
        
        #Get all objids where there are at least 3 metallicites, currently the limit has been removed
        db = sql.connect(host="ip-******.main.oberlin.edu", user="******", password="******", database="****")
        c = db.cursor()
        c.execute(f'''select distinct(objid) from dr17_metallicities where {type} is not NULL limit 6000 ;''') #putting a limit makes this way freaking faster
        galaxiesIDS = c.fetchall()
        db.close()


        totalNum = len(galaxiesIDS)
        numSoFar = 0

        print(datetime.datetime.now().time())
        print(f"Objids pulled for {type} from dr17")
        print()

        for galaxyID in galaxiesIDS:
            
            useID = galaxyID[0]

            #Gets metallicities and radii, no edits
            db = sql.connect(host="ip-******.main.oberlin.edu", user="******", password="******", database="****")
            c = db.cursor()
            c.execute(f'''select m.{type}, r.r_eff_nsa, r.agn_flag_smc from dr17_metallicities m, dr17_spaxels_uber r where m.objid='{useID}' and (m.{type} is not NULL or r.agn_flag_smc is not NULL) and m.spaxid=r.spaxID;''')
            rows = c.fetchall()
            db.close()

            #Here we are splitting the arraw into three 
            values = np.asarray(rows)
            metal = values[:,0].astype(float)
            radii = values[:, 1].astype(float)
            flags = values[:, 2]

            #STORE GALAXYID AND OTHER PIECES OF DATA IN A FILE
            universe_data[useID] = (metal, radii, flags)

            numSoFar += 1
            
            print(f"Galaxy data stored. {numSoFar} galaxies out of {totalNum} done for {type}.")
            print(f"There are {len(universe_data.keys())} keys in the dictionary rn.")
            print(datetime.datetime.now().time())
            print()


        #SAVE THE FILE
        pickle.dump(universe_data, open(f"Data/dr17 {type}.pkl", 'wb'))

    db.close()
    print(datetime.datetime.now().time())

    exit()





if __name__ == "__main__":
    main()