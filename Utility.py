import cPickle
import smtplib
from email.MIMEText import MIMEText

import scipy

eig = scipy.linalg.eig
linspace = scipy.linspace

def send_email(to_addr, from_addr, subject, message):
    """
    Send a plain-text email to a single address.
    """
    # Create a text/plain message
    msg = MIMEText(message)
    msg['Subject'] = subject
    msg['To'], msg['From'] = to_addr, from_addr

    # From the documentation for email... Not sure what it means. :-)
    ## Send the message via our own SMTP server, but don't include the
    ## envelope header.
    s = smtplib.SMTP()
    s.connect()
    s.sendmail(from_addr, [to_addr], msg.as_string())
    s.close()

def save(obj, filename):
    """
    Save an object to a file
    """
    f = file(filename, 'wb')
    cPickle.dump(obj, f, 2)
    f.close()

def load(filename):
    """
    Load an object from a file
    """
    f = file(filename, 'rb')
    obj = cPickle.load(f)
    f.close()
    return obj
