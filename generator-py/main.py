from context import SNVGeneratorContext
from model import SNVModel

if __name__ == "__main__":
    model = SNVModel(SNVGeneratorContext())
    model.generate_structure(tree_size=10)
    print("x")