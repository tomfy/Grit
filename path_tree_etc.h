// path_tree_etc.h    path_tree, Histogram, Stats stuff

void free_path_tree(Path_tree* tree);
void free_path_subtree(Path_tree_node* tree_node);
Path_tree* alloc_path_tree(void);
void print_path_tree(FILE* fp, Path_tree* tree, int print_it);
void path_tree_node_insert_multi(Path_tree_node** tree_node, State* state, int weight);
Path_tree_node*  extract_smallest_from_tree(Path_tree* tree);
Path_tree_node*  extract_smallest_from_subtree(Path_tree_node** tree_node_ptrptr);
double check_path_freq_vs_pi(Path_tree* tree);

