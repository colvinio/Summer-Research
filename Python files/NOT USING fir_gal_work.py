import pickle

def main():

    # Specify the file path of the pickle file
    file_path = "/Users/colviniorio/Desktop/Research/dr17 Z94_smc, 150 galaxies.pkl"

    # Open the pickle file in binary mode
    with open(file_path, "rb") as file:
        # Load the data from the pickle file
        data = pickle.load(file)

    objIDs = data.keys()

    print (objIDs)

if __name__ == "__main__":
    main()