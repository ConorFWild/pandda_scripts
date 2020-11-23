


class Data:
    annotation: 

class Sample:
    fragment
    polymer
    
    

class Dataset:
    fragment
    polymer
    
    def sample() -> Sample:
        ...
    
    

class Model:
    ...
    
def train():
    args: Args
    
    model: Model
    
    dataset: Dataset
    
    optimiser: Optimiser
    
    for epoch in epochs:
    
        for batch in dataset:
            annotation = batch.annotation
            input = batch.input
            
            prediction = model(input)
            
            optimiser.step()
            
        model.summary()
        
    model.save()