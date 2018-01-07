from flask import render_template
from flask import request, Markup, flash
from mmdx import app
import os
from flask import Flask, redirect, url_for
from werkzeug.utils import secure_filename
from flask import send_from_directory
import pandas as pd
import pickle
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectPercentile, f_classif
import numpy as np
import cv2
# import math
# import matplotlib.pyplot as plt
import dicom as dicom
# from skimage.transform import rescale
# from skimage import img_as_ubyte
# from io import BytesIO
# from skimage.filters import sobel
# from sqlalchemy import create_engine
# from sqlalchemy_utils import database_exists, create_database
# import pandas as pd
# import psycopg2

UPLOAD_FOLDER = '/uploads/'
RNA_EXTENSIONS = set(['txt', 'csv'])
MRI_EXTENSIONS = set(['dcm', 'tiff', 'png', 'jpg', 'jpeg'])
CLIN_EXTENSIONS = set(['txt', 'csv'])

#app = Flask(__name__)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/')
@app.route('/index')
def index():
    return render_template('jumbotron.html')


@app.route('/cool_form', methods=['GET', 'POST'])
def cool_form():
    if request.method == 'POST': # if data == [] return to jumbotron
        # do stuff when the form is submitted

        # redirect to end the POST handling
        # the redirect can be to the same route or somewhere else
        return redirect(url_for('jumbotron.html'))

    # show the form, it wasn't submitted
    return render_template('mmdxpredictearly.html'
    #, score = ENSEMBLEpredict(newdata)
    )

@app.route('/go_back', methods=['GET', 'POST'])
def go_back():
    if request.method == 'POST':
        # do stuff when the form is submitted

        # redirect to end the POST handling
        # the redirect can be to the same route or somewhere else
        return redirect(url_for('jumbotron.html'))

    # show the form, it wasn't submitted
    return render_template('jumbotron.html')

# data = RNA + MRI + clin
# newdata = transformed and normalize DATA for ensemble stacking

# RNApath = 'path/to/model/on/server'
# MRIpath = 'path/to/model/on/server'
# CLINpath = 'path/to/model/on/server'

rna_percentile_path = '/home/ubuntu/mmdxapp-master/mmdx/uploads/rna_percentile_selector.pkl'
rna_scaler_path = '/home/ubuntu/mmdxapp-master/mmdx/uploads/rna_standScaler.pkl'
rna_model_path = '/home/ubuntu/mmdxapp-master/mmdx/uploads/NB_rna_model.pkl'
rna_varthresh_path = '/home/ubuntu/mmdxapp-master/mmdx/uploads/vt_scaler.pkl'

def unpickle(filepath):
	with open(filepath, 'rb') as item:
		estimator = pickle.load(item)
	return estimator

# MAY HAVE TO DO THREE COPIES OF THIS AND THIS ONE COULD COMBINE/AVERAGE ALL 3 MODELS

# def ENSEMBLEpredict(model1, model2, model3, newdata):  
#		preprocess newdata
#		build ensemble stacking model OR just take average prediction
#		predict
#		return score

def RNApredict(rna_filepath):  
	rna = pd.read_table(rna_filepath,sep='\t', na_values='NA',low_memory=False,nrows=1 )
	rna = rna.drop(['pat.bcr','stage', 'class.y', 'class.x', 'patient.stage_event.pathologic_stage'],axis=1)
	selectov = unpickle(rna_varthresh_path) # variance threshold
	selector = unpickle(rna_percentile_path) # select percentile
	scaler = unpickle(rna_scaler_path) # standard scaler
	rna_model = unpickle(rna_model_path)
	rna = selectov.transform(rna)
	rna = selector.transform(rna)
	rna = scaler.transform(rna)
	probs = rna_model.predict_proba(rna)[:, 1]
	return rna_model, probs

# def MRIpredict(img_data):  
#		read MRI input data
#		transform input data
#		mri.model = unpickle(MRIpath)
#		predict
#		return score

# def CLINpredict(clin):  
#		newdata = read CLIN input data
#		newdata = transform input data
#		clin.model = unpickle(CLINpath)
#		predict(clin.model, newdata)
#		return score





# @app.route("/chart")
# def chart():
#     labels = ["Early","Late"]
#     values = [78,22]
#     colors = [ "#F7464A", "#46BFBD"  ]
#     return render_template('chart.html', set=zip(values, labels, colors))


# @app.route('/go')
# def go():
#     query = request.args.get('query', '')
#     return render_template(
#         'mmdxpredict.html'#,
#         query=query,
#     )


def allowed_file_rna(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in RNA_EXTENSIONS
def allowed_file_mri(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in MRI_EXTENSIONS
def allowed_file_clin(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in CLIN_EXTENSIONS

@app.route('/rna_input', methods=['GET', 'POST'])
def rna_input():
	imgID = request.form.get('rna_input')
		# print(str(imgID))

		# Match subject ID to image file names stored on the data base
	if str(imgID) == 'Patient1_RNASeq':
		filename = 'RNAseq2.txt'
	elif str(imgID) == 'Patient2_RNASeq':
		filename = 'RNAseq2.txt'
	elif str(imgID) == 'Patient3_RNASeq':
		filename = 'RNAseq2.txt'
	else:
		filename = 'RNAseq2.txt'
	fullpath = '/home/ubuntu/mmdxapp-master/mmdx/uploads/' + filename
	#data = pd.read_csv(fullpath)
	rna_model, probs = RNApredict(fullpath)
	return render_template('jumbotron_rna_done.html', score = probs)
#     if request.method == 'POST':
#         # check if the post request has the file part
#         if 'file' not in request.files:
#             flash('No file part')
#             return redirect(request.url)
#         file = request.files['file']
#         # if user does not select file, browser also
#         # submit a empty part without filename
#         if file.filename == '':
#             flash('No selected file')
#             return redirect(request.url)
#         if file and allowed_file_rna(file.filename):
#             filename = secure_filename(file.filename)
#             filename = 'RNAseq2.txt'
#             file.save('/home/ubuntu/mmdxapp-master/mmdx/uploads/' + filename)
#             #rna_path = '/home/ubuntu/mmdxapp-master/mmdx/uploads/' + filename
#             #probs = RNApredict(rna_path)
#             #file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
#             return render_template('jumbotron_rna_done.html'
#             , score = 'File Upload Successful')
    # return '''
#     <!doctype html>
#     <title>Upload new File</title>
#     <h1>Upload new File</h1>
#     <form method=post enctype=multipart/form-data>
#       <p><input type=file name=file>
#          <input type=submit value=Upload>
#     </form>
#     '''


@app.route('/mri_input', methods=['GET', 'POST'])
def mri_input():
	imgID = request.form.get('rna_input')

	if str(imgID) == 'Patient 1 MRI':
		filename = 'dicom_dir'
	elif str(imgID) == 'Patient 2 MRI':
		filename = 'dicom_dir'
	elif str(imgID) == 'Patient 2 MRI':
		filename = 'dicom_dir'
	else:
		filename = 'dicom_dir'
	fullpath = '/home/ubuntu/mmdxapp-master/mmdx/uploads/' + filename
	files = os.listdir(fullpath)
	scan = dicom.read_file(fullpath + '/' + files[0])
	
	IMG_PX_SIZE= 150
	new_scan = cv2.resize(np.array(scan.pixel_array),(IMG_PX_SIZE,IMG_PX_SIZE))
	scan = np.array(scan.pixel_array)
	#scan /= 255.
	
	#score = scan.shape
	mean = np.mean(new_scan)
	
	# Plot out the image and save it as a png file
# 	plt.imshow(scan, origin='lower', cmap=plt.cm.gray)
# 	plt.axis('off')
# 	figfile = BytesIO()
# 	plt.savefig(figfile, format='png')
# 	figfile.seek(0)    # Rewind to beginning of file
# 	figdata_png = figfile.getvalue()    # Extract string

	# Convert the data to base64 format
	#import base64
	#figdata_png = base64.b64encode(figdata_png).decode('ascii')
	return render_template('jumbotron_mri_done.html', score = mean)
    
#     if request.method == 'POST':
#         # check if the post request has the file part
#         if 'file' not in request.files:
#             flash('No file part')
#             return redirect(request.url)
#         file = request.files['file']
#         # if user does not select file, browser also
#         # submit a empty part without filename
#         if file.filename == '':
#             flash('No selected file')
#             return redirect(request.url)
#         if file and allowed_file_mri(file.filename):
#             filename = secure_filename(file.filename)
#             filename = 'RNAseq2.txt'
#             file.save('/home/ubuntu/mmdxapp-master/mmdx/uploads/' + filename)
#             #file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
#             return render_template('jumbotron_mri_done.html'
#             , score = 'File Upload Successful')
#             #return 'successful'
#             #return redirect(url_for('uploaded_file'))
#             #return redirect(url_for('uploaded_file',filename=filename))
#     return '''
#     <!doctype html>
#     <title>Upload new File</title>
#     <h1>Upload new File</h1>
#     <form method=post enctype=multipart/form-data>
#       <p><input type=file name=file>
#          <input type=submit value=Upload>
#     </form>
#     '''

@app.route('/clin_input', methods=['GET', 'POST'])
def clin_input():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit a empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file_clin(file.filename):
            filename = secure_filename(file.filename)
            filename = 'RNAseq2.txt'
            file.save('/home/ubuntu/mmdxapp-master/mmdx/uploads/' + filename)
            #file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            return render_template('jumbotron_clin_done.html'
            , score = 'File Upload Successful')
            #return 'successful'
            #return redirect(url_for('uploaded_file'))
            #return redirect(url_for('uploaded_file',filename=filename))
    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Upload>
    </form>
    '''


    
#@app.route('/uploads/<filename>')
#def uploaded_file(filename):
#	return 'successful'
	#return render_template('mmdxpredict.html')
    #return send_from_directory(app.config['UPLOAD_FOLDER'],filename)