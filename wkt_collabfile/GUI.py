import tkinter as tk
from tkinter import filedialog, messagebox

class GlycanAnalysisApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Glycan Analysis Application")
        
        # Page 1: Home Page
        self.home_frame = tk.Frame(self.root)
        self.home_frame.pack()
        
        tk.Label(self.home_frame, text="Select to Train Your Profile or Predict Your Sample").pack(pady=10)
        self.train_button = tk.Button(self.home_frame, text="Train", command=self.train_mode)
        self.train_button.pack(pady=5)
        
        self.predict_button = tk.Button(self.home_frame, text="Predict", command=self.predict_mode)
        self.predict_button.pack(pady=5)
        
    def clear_frame(self):
        for widget in self.root.winfo_children():
            widget.destroy()
    
    def train_mode(self):
        # Page 2: File Upload for Training
        self.clear_frame()
        self.train_file_path = None
        
        tk.Label(self.root, text="File Upload Page - Train Mode").pack(pady=10)
        self.upload_button = tk.Button(self.root, text="Upload Training Data File", command=self.upload_train_file)
        self.upload_button.pack(pady=5)
    
    def upload_train_file(self):
        self.train_file_path = filedialog.askopenfilename()
        if self.train_file_path:
            messagebox.showinfo("File Upload", "Training Data File Uploaded Successfully")
            self.train_parameters()
    
    def train_parameters(self):
        # Page 3: Train Parameter Settings
        self.clear_frame()
        
        tk.Label(self.root, text="Train Parameter Settings Page").pack(pady=10)
        
        # Column Drop
        tk.Label(self.root, text="Column Drop").pack(pady=5)
        self.column_drop_vars = {"RT": tk.IntVar(), "filename": tk.IntVar(), "GlycoPost_ID": tk.IntVar()}
        for column in self.column_drop_vars:
            tk.Checkbutton(self.root, text=column, variable=self.column_drop_vars[column]).pack(anchor='w')
        
        # Mini Class Element
        tk.Label(self.root, text="Mini Class Element (Minimum Sample Size) - Default 5").pack(pady=5)
        self.mini_class_element = tk.Entry(self.root)
        self.mini_class_element.insert(0, "5")
        self.mini_class_element.pack(pady=5)
        
        # Test Size
        tk.Label(self.root, text="Test Size (Training/Test Split)").pack(pady=5)
        self.test_size = tk.Scale(self.root, from_=60, to=90, orient=tk.HORIZONTAL)
        self.test_size.set(80)
        self.test_size.pack(pady=5)
        
        # Cross-validation Fold
        tk.Label(self.root, text="Cross-validation Fold").pack(pady=5)
        self.cv_fold = tk.Scale(self.root, from_=3, to=10, orient=tk.HORIZONTAL)
        self.cv_fold.set(5)
        self.cv_fold.pack(pady=5)
        
        # Model Select
        tk.Label(self.root, text="Model Select").pack(pady=5)
        self.models = {"SVM": tk.IntVar(), "RandForest": tk.IntVar(), "KNN": tk.IntVar(), "XGBoost": tk.IntVar()}
        for model in self.models:
            tk.Checkbutton(self.root, text=model, variable=self.models[model]).pack(anchor='w')
        
        # Submit Button
        self.submit_button = tk.Button(self.root, text="Submit", command=self.submit_training)
        self.submit_button.pack(pady=10)
    
    def submit_training(self):
        # Collect training parameters and proceed (not implemented)
        messagebox.showinfo("Training", "Training Started with Selected Parameters")
    
    def predict_mode(self):
        # Page 2: File Upload for Prediction
        self.clear_frame()
        self.predict_file_path = None
        
        tk.Label(self.root, text="File Upload Page - Predict Mode").pack(pady=10)
        self.upload_button = tk.Button(self.root, text="Upload Prediction Data File", command=self.upload_predict_file)
        self.upload_button.pack(pady=5)
    
    def upload_predict_file(self):
        self.predict_file_path = filedialog.askopenfilename()
        if self.predict_file_path:
            messagebox.showinfo("File Upload", "Prediction Data File Uploaded Successfully")
            self.predict_parameters()
    
    def predict_parameters(self):
        # Page 3: Predict Parameter Settings
        self.clear_frame()
        
        tk.Label(self.root, text="Predict Parameter Settings Page").pack(pady=10)
        
        # Column Drop
        tk.Label(self.root, text="Column Drop").pack(pady=5)
        self.column_drop_vars = {"RT": tk.IntVar(), "filename": tk.IntVar(), "GlycoPost_ID": tk.IntVar()}
        for column in self.column_drop_vars:
            tk.Checkbutton(self.root, text=column, variable=self.column_drop_vars[column]).pack(anchor='w')
        
        # Submit Button
        self.submit_button = tk.Button(self.root, text="Submit", command=self.submit_prediction)
        self.submit_button.pack(pady=10)
    
    def submit_prediction(self):
        # Collect prediction parameters and proceed (not implemented)
        messagebox.showinfo("Prediction", "Prediction Started with Selected Parameters")

if __name__ == "__main__":
    root = tk.Tk()
    app = GlycanAnalysisApp(root)
    root.mainloop()