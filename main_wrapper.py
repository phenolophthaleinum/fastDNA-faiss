import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="fastDNA-faiss")
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory with virus samples of vectors.")
    parser.add_argument("-o", "--output", required=True,
                        help="Path to result file with rankings.")
    parser.add_argument("-d", "--dim", required=True,
                        help="Dimensionality of vectors")
    parser.add_argument("--n_samples", required=True,
                        help="Number of samples from each host genome")
    parser.add_argument("-k", "--k_nearest", required=True,
                        help="How many (k) nearest samples should be found")
    parser.add_argument("-f", "--faiss_index", required=True,
                        help="Path to faiss index")
    parser.add_argument("-m", "--map", required=True,
                        help="Path to sample map")

    args = parser.parse_args()

    global_rank = {}
    run_procedure(args.input_dir, args.output, int(args.k_nearest), int(args.dim), int(args.n_samples),
                  args.faiss_index, args.map)