#include <iostream>
#include <unordered_map>
#include <fstream>
#include <trees/PLLUnrootedTree.hpp>


struct UnorderedPairInt {
  UnorderedPairInt():first(0), second(0) {}
  UnorderedPairInt(unsigned int f, unsigned int s): first(f), second(s) {}
  bool operator ==(const UnorderedPairInt &other) const {
    return (first == other.first && second == other.second) ||
      (first == other.second && second == other.first);
  }
  unsigned int first;
  unsigned int second;
};

namespace std {
  template<> struct hash<UnorderedPairInt> {
    size_t operator()(const UnorderedPairInt &p) const {
      return p.first ^ p.second;
    }
  };
}

class CompressedTreeList {
public:
  CompressedTreeList(const std::string &treeListFile) {
    std::ifstream is(treeListFile);
    std::string line;
    unsigned int treeCount = 0;
    while (std::getline(is, line)) {
      if (line.size() < 2) {
        continue;
      }
      PLLUnrootedTree tree(line, false);
      treeCount++;
      if (_labelToId.size() == 0) {
        // first tree: fill the leaf<->ID mapping
        for (auto &label: tree.getLeavesLabels()) {
          _idToLabel.push_back(label);
          _labelToId.insert({label, _labelToId.size()});
          UnorderedPairInt p;
          // add a fake subtree for leaves
          _subtreeChildren.push_back(p); 
        }
        _maxID = _labelToId.size() - 1;
      }
      // arbitrarily root all the tree at the same leaf branch
      auto root = tree.findLeaf(_idToLabel[0])->back;
      auto treeID = mapSubtree(root);
      if (_weightedTrees.find(treeID) == _weightedTrees.end()) {
        _weightedTrees.insert({treeID, 0});
      }
      _weightedTrees[treeID]++;
    }
    std::cout << "Number of input trees " << treeCount << std::endl;
    std::cout << "Number of unique trees " << _weightedTrees.size() << std::endl;
    std::cout << "Number of leaves: " << _labelToId.size() << std::endl;
    std::cout << "Number of different subtrees: " << _subtrees.size() << std::endl;
  }

  void dumpUncompressed(const std::string &output) {
    std::ofstream os(output);
    for (auto &weightedTree: _weightedTrees) {
      auto treeID = weightedTree.first;
      auto treeCount = weightedTree.second;
      auto newick = getNewick(treeID);
      for (unsigned int i = 0; i < treeCount; ++i) {
        os << newick << std::endl;
      }
    }
  }

  void dumpCompressed(const std::string &output) {
    std::ofstream os(output, std::ios::binary);
    // write magic string
    std::string magic = "magicctl";
    os.write(magic.c_str(), magic.size());
    // write labels
    unsigned int leafNumber = _idToLabel.size();
    writeUInt(leafNumber, os);
    for (auto &label: _idToLabel) {
      writeString(label, os);
    }
    // write subtree identifiers
    writeUInt(_subtreeChildren.size(), os);
    for (unsigned int i = 0; i < _subtreeChildren.size(); ++i) {
      writeUInt(i, os);
      writeUInt(_subtreeChildren[i].first, os);
      writeUInt(_subtreeChildren[i].second, os);
    }
  }

  void loadCompressed(const std::string &inputFile) {
    std::ifstream is(inputFile, std::ios::binary);
  }

private:
  static void writeUInt(unsigned int i, std::ofstream &os) {
    os.write((char *)&i, sizeof(unsigned int));
  }

  static void writeString(const std::string &str, std::ofstream &os) {
    writeUInt(str.size(), os);
    os.write(str.c_str(), str.size());
  }
  
  // todo: we could cache the aux newicks
  std::string getNewickAux(unsigned int id) {
    if (id < _idToLabel.size()) {
      // leaf
      return _idToLabel[id];
    } 
    // internal node
    std::string newick = "(";
    auto p = _subtreeChildren[id];
    newick += getNewickAux(p.first);
    newick += ",";
    newick += getNewickAux(p.second);
    newick += ")";
    return newick;
  };

  std::string getNewick(unsigned int treeID) {
    auto p = _subtreeChildren[treeID];
    std::string newick = "(";
    // add the (arbitrary) root label
    newick += _idToLabel[0]; 
    newick += ",";
    newick += getNewickAux(p.first);
    newick += ",";
    newick += getNewickAux(p.second);
    newick += ");";
    return newick;
  }


  /**
   * map the subtree to an ID and return this ID
   */
  unsigned int mapSubtree(pll_unode_t *node) {
    if (!node->next) { // terminal node
      return _labelToId[node->label];
    }
    auto child1 = node->next->back;
    auto child2 = node->next->next->back;
    auto id1 = mapSubtree(child1);
    auto id2 = mapSubtree(child2);
    return addSubtreeAndGetIndex(id1, id2);
  }
  
  /**
   * aux function for mapSubtree
   */
  unsigned int addSubtreeAndGetIndex(unsigned int id1, unsigned int id2) {
    UnorderedPairInt p(id1, id2);
    auto it = _subtrees.find(p);
    if (it == _subtrees.end()) {
      // new subtree, add it
      _maxID++;
      _subtrees.insert({p, _maxID});
      _subtreeChildren.push_back(p);
      return _maxID;
    } else {
      // this subtree already exists
      return it->second;
    }
  }

  std::unordered_map<std::string, unsigned int> _labelToId;
  std::vector <std::string> _idToLabel; 
  std::unordered_map<UnorderedPairInt, unsigned int> _subtrees;
  std::vector<UnorderedPairInt> _subtreeChildren;
  std::unordered_map<unsigned int, unsigned int> _weightedTrees;
  unsigned int _maxID;
};

int main(int argc, char** argv)
{
  if (argc != 3) {
    std::cerr << "Syntax: tree_list_compress tree_list_file output" << std::endl;
    return 1;
  }
  std::string treeListFile(argv[1]);  
  std::string output(argv[2]);
  CompressedTreeList treeList(treeListFile);
  treeList.dumpUncompressed("plop1nobl.newick");
  treeList.dumpCompressed(output);

  return 0;
}


