# Week 9 - Copier project templates

This week we will look at using a project template, Copier, to turn the CDK-Depict code you wrote last week into a simple Python package.

**What is Copier?**
   - Copier is a utility tool created to facilitate the duplication and customization of project templates. It is written in Python and enables you to efficiently generate projects by copying a template and making specified alterations defined in an easy configuration file (often YAML).
   - This tool is designed to help reduce redundant setup tasks and ensures consistency across multiple projects by allowing users to create a single project template that can be reused and customized without starting from scratch each time.

Lets go through the steps required to make this package.

#### Step 1: Set Up Your Environment

```bash
conda activate ppchem
```

#### Step 2: Organise the code

Go through the Jupyter notebook from last week and ensure that your code is neatly organised into functions. Once this is done move the code into a new nodebook which should contain only the code required to query the CDK-Depict website. 
Place this notebook in a new folder called 'CDK-Package'.

Remember you can use the following commands to create folders and move files between directories. For the *mv* command to work you must navigate your terminal to the directory in which the file you want to move is located. In your case it should be something like '/path/to/your/practical-programming-in-chemistry-exercises/week_08/your_new_notebook.ipynb'

```bash
mkdir CDK-Package
cd /path/to/your/practical-programming-in-chemistry-exercises/week_08/your_new_notebook.ipynb
mv your_new_notebook.ipynb /path/to/your/CDK-Package
```

#### Step 3: Convert Your Notebook

While Jupyter notebooks are excellent tools for quickly developing small snippets of code, they are not suited for use in a Python package. To make your code compatible with a Python package we
Export your Jupyter Notebook (`YourNotebook.ipynb`) to a standard Python script (`module.py`). You can do this directly from the Jupyter interface or use `nbconvert`:

```bash
jupyter nbconvert --to script your_new_notebook.ipynb
```

This creates `YourNotebook.py`. Rename and organize this file into a meaningful module name.

#### Step 4: Using Copier to Create the Package Structure**

Instead of manually creating directories and files as outlined previously, you can use the Copier template designed for consistent and quick setup of scientific projects/folders. Here’s how you can do it:

1. **Install Copier**:
   First, ensure Copier is installed in your environment:

   ```bash
   pip install copier
   ```

2. **Generate Project Structure**:
   Utilize Copier with the schwallergroup template to create your project structure. Make sure you're still within your `ppchem` environment:

   ```bash
   copier https://github.com/schwallergroup/copier-liac.git /path/to/you/CDK-Package
   ```

   Enter the directory path where you want your new project to be initiated. Follow the on-screen prompts provided by Copier to customize the project (such as naming modules or defining author details).

3. **Place Your Converted Module**:
   Move or copy the Python scripts (converted from your Jupyter notebooks in Step 3) into the appropriate directories within this newly created package structure. Typically, this will be under the main package directory.

Proceed with the next steps as previously outlined:

#### Step 6 (updated): Initialize Git Repository
Now, initialize a git repository to start version control within the newly created directory structure:

```bash
cd /path/to/you/CDK-Package
git init
git add .
git commit -m "Initial package setup with Copier"
```

#### Steps 7 to 9
Follow the original steps to install your package. These actions involve using `pip install -e .`.

#### Advantages of Using Copier:
- **Standardized Setup**: Ensures all projects start with a consistent, error-free base.
- **Efficiency**: Saves time and reduces manual errors in setting up a package structure.
- **Scalability**: Easy to update and scale your project structure as Copier templates can be modified and reused.

### Conclusion
By incorporating Copier in the setup process, you enhance reproducibility and consistency across scientific programming projects in chemistry. This approach not only streamlines the workflow but also guarantees that all necessary best practices in software development are followed from the outset.

