##### MAIN ####
import load_data
import heuristics_PSSM
import hobohm1_centroids
import get_hobohm_PSSM
import predictions


load_data.main()
heuristics_PSSM.main()
get_hobohm_PSSM.main()
#hobohm1_centroids.main()
predictions.main()


print("Done")

