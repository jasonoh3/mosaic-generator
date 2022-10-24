/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    if (first[curDim] < second[curDim]) {
      return true;
    } else if (first[curDim] > second[curDim]) {
      return false;
    } 
    
    if (first < second) {
      return true;
    }
    return false;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    double PTDist = 0.0;
    double CBTDist = 0.0;

    for (int i = 0; i < Dim; ++i) {
      PTDist += (target[i] - potential[i]) * (target[i] - potential[i]);
      CBTDist += (target[i] - currentBest[i]) * (target[i] - currentBest[i]);
    }

    if (PTDist < CBTDist) {
      return true;
    } else if (PTDist > CBTDist) {
      return false;
    }

    if(potential < currentBest) {
      return true;
    }
    return false;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    size = 0;

    if (newPoints.empty()) {
      root = NULL;
    } else {
      vector<Point<Dim>> copiedPoints = newPoints;
      buildTree(copiedPoints, 0, 0, copiedPoints.size() - 1, root);
    }
}

template <int Dim>
int KDTree<Dim>::partition(vector<Point<Dim>>& list, int dim, int left, int right, int pivotIndex) {
  Point<Dim> pivotValue = list[pivotIndex];
  list[pivotIndex] = list[right];
  list[right] = pivotValue;

  int storeIndex = left;
  for (int i = left; i < right; ++i) {
    if (smallerDimVal(list[i], pivotValue, dim)) {
      Point<Dim> temp = list[storeIndex];
      list[storeIndex] = list[i];
      list[i] = temp;
      ++storeIndex;
    }
  }

  Point<Dim> temp = list[right];
  list[right] = list[storeIndex];
  list[storeIndex] = temp;

  return storeIndex;
}

template <int Dim>
Point<Dim> KDTree<Dim>::quickSelect(vector<Point<Dim>>& list, int dim, int left, int right, int k) {
  for (unsigned i = 0; i < list.size(); ++i) {
    if (left == right) {
      return list[left];
    }
    
    int pivotIndex = k;
    pivotIndex = partition(list, dim, left, right, pivotIndex);

    if (k == pivotIndex) {
      return list[k];
    } else if (k < pivotIndex) {
      right = pivotIndex - 1;
    } else {
      left = pivotIndex + 1;
    }
  }
  
  return NULL;
}

template <int Dim>
void KDTree<Dim>::buildTree(vector<Point<Dim>>& points, int dim, int left, int right, KDTreeNode*& curRoot) {
  if (left <= right) {
    int middle = (left + right) / 2;
    curRoot = new KDTreeNode(quickSelect(points, dim, left, right, middle));
    ++size;

    buildTree(points, (dim + 1) % Dim, left, middle - 1, curRoot->left);
    buildTree(points, (dim + 1) % Dim, middle + 1, right, curRoot->right);
  }
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  copy(this->root, other.root);
  size = other.size;
}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode* current, KDTreeNode* other) {
  if (other == NULL) {
    return;
  }
  current = new KDTreeNode(other->point);
  copy(current->left, other->left);
  copy(current->right, other->right);
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  if (&rhs != this) {
    destroy(root);
    copy(this->root, rhs.root);
    size = rhs.size;
  }

  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  destroy(root);
}

template <int Dim>
void KDTree<Dim>::destroy(KDTreeNode*& root) {
  if (root == NULL) {
    return;
  }
  destroy(root->left);
  destroy(root->right);
  delete root;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const {
  return findNearestNeighbor(query, 0, root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query, int dim, const KDTreeNode* curRoot) const {
  if (curRoot->left == NULL && curRoot->right == NULL) {
    return curRoot->point;
  }

  Point<Dim> nearest;
  bool wentRight = true;
  if (smallerDimVal(query, curRoot->point, dim) && curRoot->left != NULL) {
    nearest = findNearestNeighbor(query, (dim + 1) % Dim, curRoot->left);
    wentRight = false;
  } else {
    nearest = findNearestNeighbor(query, (dim + 1) % Dim, curRoot->right);
  }

  if (shouldReplace(query, nearest, curRoot->point)) {
    nearest = curRoot->point;
  }

  double radius = 0.0;
  for (int i = 0; i < Dim; ++i) {
    radius += (nearest[i] - query[i]) * (nearest[i] - query[i]);
  }
  
  double splitDist = (curRoot->point[dim] - query[dim]) * (curRoot->point[dim] - query[dim]);
  

  Point<Dim> tempNearest;
  if (radius >= splitDist) {
    if (wentRight && curRoot->left != NULL) {
      tempNearest = findNearestNeighbor(query, (dim + 1) % Dim, curRoot->left);
    } else {
      tempNearest = findNearestNeighbor(query, (dim + 1) % Dim, curRoot->right);
    }
    if (shouldReplace(query, nearest, tempNearest)) {
      nearest = tempNearest;
    }
  }

  return nearest;
}

